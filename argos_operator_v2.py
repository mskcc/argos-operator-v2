import json
import logging
import os
import re
from datetime import datetime as dt

from voyager_sdk.file_repository import FileRepository
from voyager_sdk.operator.operator import Operator
from voyager_sdk.protocols.processors.file_processor import FileProcessor

LOGGER = logging.getLogger(__name__)


PDX_SPECIMEN_TYPES = ["pdx", "xenograft", "xenograftderivedcellline"]
NON_PDX_SPECIMEN_TYPES = [
    "biopsy",
    "blood",
    "cellLine",
    "cfdna",
    "fingernails",
    "nonpdx",
    "normal",
    "organoid",
    "other",
    "rapidautopsy",
    "resection",
    "saliva",
    "tumor",
    "poolednormal",
    "dmp",
]


class ArgosOperatorV2(Operator):
    """Base operator class for the SDK"""

    def __init__(
        self,
        request_id=None,
        runs=[],
        pipeline=None,
        pairing=None,
        output_directory_prefix=None,
        file_group=None,
        job_group_id=None,
        job_group_notifier_id=None,
        **kwargs,
    ):
        """
        request_id: metadata key:igoRequestId
        runs: runs[]
        pipeline: {
            "pipeline_format": "",
            "pipeline_link": "",
            "pipeline_version: "",
            "pipeline_entrypoint": ""
        },
        file_group: file_group_id
        """
        super().__init__(
            request_id,
            runs,
            pipeline,
            pairing,
            output_directory_prefix,
            file_group,
            job_group_id,
            job_group_notifier_id,
        )

    def get_samples_from_data(self, data):
        samples = list()
        # group by igoId
        igo_id_group = dict()
        for sample in data:
            igo_id = sample["metadata"]["primaryId"]
            if igo_id not in igo_id_group:
                igo_id_group[igo_id] = list()
            igo_id_group[igo_id].append(sample)

        for igo_id in igo_id_group:
            sample = igo_id_group[igo_id][0]
            sample_name = sample["metadata"]["cmoSampleName"]
            if "poolednormal" in sample_name.lower():
                samples.append(
                    self.build_sample(
                        igo_id_group[igo_id], ignore_sample_formatting=True
                    )
                )
            else:
                samples.append(self.build_sample(igo_id_group[igo_id]))
        return samples

    def build_data_list(self, files):
        data = list()
        for f in files:
            sample = dict()
            sample["id"] = f.file_id
            sample["path"] = f.path
            sample["file_name"] = f.file_name
            sample["metadata"] = f.metadata
            data.append(sample)
        return data

    def get_jobs(self):
        argos_jobs = list()
        dmp_samples = list()
        files = self.get_files(self.request_id)

        data = self.build_data_list(files)

        samples = self.get_samples_from_data(data)
        argos_inputs, error_samples = self.construct_argos_jobs(
            samples, logger=self.logger
        )
        sample_pairing = self.get_pairing_from_argos_inputs(argos_inputs)
        sample_mapping, filepaths = self.get_mapping_from_argos_inputs(argos_inputs)
        argos_jobs = self.get_argos_jobs(argos_inputs)

        return argos_jobs

    def remove_with_caveats(self, samples):
        """
        Removes samples from a list of samples if they either
        don't contain a 'sampleNameMalformed', which happens when function
        format_sample_name returns it
        """
        data = list()
        error_data = list()
        for sample in samples:
            add = True
            sample_id = sample["sample_id"]
            sample_name = sample["SM"]
            patient_id = sample["patient_id"]
            if sample_name == "sampleNameMalformed":
                add = False
                LOGGER.info(
                    "Sample name is malformed for for %s; removing from set", sample_id
                )
            if not patient_id:
                add = False
                LOGGER.info("No patient ID for sample %s; removing from set", sample_id)
            elif isinstance(patient_id, str):
                if not patient_id.strip():
                    add = False
                    LOGGER.info(
                        "Empty string for patient ID in sample %s; removing from set",
                        sample_id,
                    )
            if add:
                data.append(sample)
            else:
                error_data.append(sample)
        return data, error_data

    def construct_argos_jobs(self, samples, logger=None):
        samples, error_samples = self.remove_with_caveats(samples)
        pairs = self.compile_pairs(samples)
        number_of_tumors = len(pairs["tumor"])
        argos_jobs = list()
        for i in range(0, number_of_tumors):
            tumor = pairs["tumor"][i]
            normal = pairs["normal"][i]
            project_id = tumor["request_id"]
            assay = tumor["bait_set"]
            patient_id = tumor["patient_id"]
            pi = tumor["pi"]
            pi_email = tumor["pi_email"]
            job = dict()
            tumor_specimen_type = normalize_igo_text_field(
                pairs["tumor"][i]["specimen_type"]
            )
            normal_sample = format_sample(normal, tumor_specimen_type)
            tumor_sample = format_sample(tumor, tumor_specimen_type)
            job["tumor"] = tumor_sample
            job["normal"] = normal_sample
            job["assay"] = assay
            job["pi"] = pi
            job["pi_email"] = pi_email
            job["patient_id"] = patient_id
            references = convert_references(project_id, assay, pi, pi_email)
            job.update(references)
            argos_jobs.append(job)
        return argos_jobs, error_samples

    def get_files(self, request_id):
        files = FileRepository.filter(metadata=[f"igoRequestId:{request_id}"])
        # cnt_tumors = FileRepository.filter(
        #     metadata=[
        #         f"igoRequestId:{request_id}"
        #         "tumorOrNormal:Tumor"
        #     ]
        # ).count
        return files

    def build_sample(self, data, ignore_sample_formatting=False):
        """
        Given some data - which is a list of samples originally from the LIMS, split up into one file
        per index - the data is then compiled into one sample dictionary consisting of one or more
        pairs of fastqs

        Note that ID and SM are different field values in ARGOS (RG_ID and ID, respectively, in ARGOS)
        but standardizing it here with what GATK sets bam headers to
        """

        samples = dict()

        for value in data:
            fpath = value["path"]
            curr_file = get_file(fpath)
            meta = value["metadata"]
            bid = value["id"]
            sequencing_center = meta["sequencingCenter"]
            platform = meta["platform"]
            request_id = meta["igoRequestId"]
            sample_id = meta["primaryId"]
            library_id = meta["libraryIgoId"]
            bait_set = meta["baitSet"]
            tumor_type = meta["tumorOrNormal"]
            specimen_type = meta["sampleClass"]
            species = meta["species"]
            cmo_sample_name = format_sample_name(
                meta["cmoSampleName"], specimen_type, ignore_sample_formatting
            )
            if cmo_sample_name == "sampleNameMalformed":
                LOGGER.error("sampleName for %s is malformed", sample_id)
            flowcell_id = meta["flowCellId"]
            barcode_index = meta["barcodeIndex"]
            cmo_patient_id = meta["cmoPatientId"]
            platform_unit = flowcell_id
            run_date = meta["runDate"]
            r_orientation = meta["R"]
            pi_name = meta.get("pi", meta["labHeadName"])
            pi_email = meta.get("pi_email", meta["labHeadEmail"])
            run_id = meta["runId"]
            preservation_type = meta["preservation"]
            rg_id = cmo_sample_name + "_1"
            run_mode = get_run_mode(meta["runMode"])
            if barcode_index:
                platform_unit = "_".join([flowcell_id, barcode_index])
            try:
                rg_id = "_".join([cmo_sample_name, platform_unit])
            except:
                LOGGER.info("RG ID couldn't be set.")
                LOGGER.info("Sample ID %s; patient ID %s", sample_id, cmo_patient_id)
                LOGGER.info(
                    "SampleName %s; platform unit %s", cmo_sample_name, platform_unit
                )
            if sample_id not in samples:
                samples[sample_id] = dict()
                sample = dict()
                sample["CN"] = sequencing_center
                sample["PL"] = platform
                sample["PU"] = list()
                sample["LB"] = library_id
                sample["tumor_type"] = tumor_type
                sample["SM"] = cmo_sample_name
                sample["species"] = species
                sample["patient_id"] = cmo_patient_id
                sample["bait_set"] = bait_set
                sample["sample_id"] = sample_id
                sample["run_date"] = run_date
                sample["specimen_type"] = specimen_type
                sample["request_id"] = request_id
                sample["pi"] = pi_name
                sample["pi_email"] = pi_email
                sample["run_id"] = run_id
                sample["preservation_type"] = preservation_type
                sample["ID"] = list()
                sample["R1"] = list()
                sample["R1_bid"] = list()
                sample["R2"] = list()
                sample["R2_bid"] = list()
                sample["fastqs"] = list()
                sample["run_mode"] = run_mode
            else:
                sample = samples[sample_id]

            # Queueing up fastqs for pairing later; RG ID and PU
            # will be assigned based on Fastqs object
            if "R1" in r_orientation or "R2" in r_orientation:
                sample["fastqs"].append(curr_file)
            else:
                # DMP bams found; assigning RG ID and PU here
                # There will always be only one DMP bam, so assign explicitly
                sample["bam"] = fpath
                sample["bam_bid"] = bid
                sample["PU"] = platform_unit
                sample["ID"] = rg_id
            samples[sample_id] = sample

        result = dict()
        result["CN"] = list()
        result["PL"] = list()
        result["PU"] = list()
        result["LB"] = list()
        result["tumor_type"] = list()
        result["ID"] = list()
        result["SM"] = list()
        result["species"] = list()
        result["patient_id"] = list()
        result["bait_set"] = list()
        result["sample_id"] = list()
        result["run_date"] = list()
        result["specimen_type"] = list()
        result["R1"] = list()
        result["R2"] = list()
        result["R1_bid"] = list()
        result["R2_bid"] = list()
        result["bam"] = list()
        result["bam_bid"] = list()
        result["request_id"] = list()
        result["pi"] = list()
        result["pi_email"] = list()
        result["run_id"] = list()
        result["preservation_type"] = list()
        result["run_mode"] = list()

        for sample_id in samples:
            sample = samples[sample_id]
            for key in sample:
                if key == "fastqs":
                    if sample["fastqs"]:
                        fastqs = Fastqs(sample["SM"], sample["fastqs"])
                        result["R1"] = fastqs.r1
                        result["R1_bid"] = fastqs.r1_bids
                        result["R2"] = fastqs.r2
                        result["R2_bid"] = fastqs.r2_bids
                        result["PU"] = fastqs.pu
                        result["ID"] = fastqs.rg_id
                else:
                    result[key].append(sample[key])
        result = self.check_and_return_single_values(result)
        return result

    def get_argos_jobs(self, argos_inputs):
        argos_jobs = list()
        number_of_inputs = len(argos_inputs)
        for i, job in enumerate(argos_inputs):
            tumor_sample_name = job["tumor"]["ID"]
            normal_sample_name = job["normal"]["ID"]

            name = "ARGOS %s, %i of %i" % (self.request_id, i + 1, number_of_inputs)
            assay = job["assay"]
            pi = job["pi"]
            pi_email = job["pi_email"]

            tags = {
                "igoRequestId": self.request_id,
                "sampleNameTumor": tumor_sample_name,
                "sampleNameNormal": normal_sample_name,
                "labHeadName": pi,
                "labHeadEmail": pi_email,
            }
            # pipeline = self.get_pipeline_id()
            log_prefix = f"{tumor_sample_name}_{normal_sample_name}"
            # log_directory = self.get_log_directory()
            if self.output_directory_prefix:
                tags["output_directory_prefix"] = self.output_directory_prefix
            argos_jobs.append(
                dict(
                    app=self._pipeline,
                    inputs=job,
                    name=name,
                    tags=tags,
                    log_prefix=log_prefix,
                )
            )
        return argos_jobs

    def check_and_return_single_values(self, data):
        """
        data is a dictionary; each key contains a list of values.

        single_values are the expected keys that should contain only one value

        Concatenating pi and pi_email AND formatting the LB field are workarounds
        because some samples would have multiple values for these but the sample dict
        it returns must have one value only in order for the pipeline to execute
        """
        single_values = [
            "CN",
            "PL",
            "SM",
            "bait_set",
            "patient_id",
            "species",
            "tumor_type",
            "sample_id",
            "specimen_type",
            "request_id",
            "run_mode",
        ]

        for key in single_values:
            value = set(data[key])
            if len(value) == 1:
                data[key] = value.pop()
            else:
                LOGGER.error("Expected only one value for %s!", key)
                LOGGER.error("Check import, something went wrong.")

        # concatenating pi and pi_email
        if data["pi"] == [None]:
            data["pi"] = [""]
        if data["pi_email"] == [None]:
            data["pi_email"] = [""]
        data["pi"] = "; ".join(set(data["pi"]))
        data["pi_email"] = "; ".join(set(data["pi_email"]))

        # hack; formats LB field so that it is a string
        library_id = [i for i in data["LB"] if i]
        number_of_library_ids = len(library_id)
        if number_of_library_ids > 0:
            library_id.sort()
            data["LB"] = "_and_".join(library_id)
        else:
            data["LB"] = data["SM"] + "_1"

        # run_ids need to be one list
        run_ids = data["run_id"]
        if any(isinstance(run_id, list) for run_id in run_ids):
            run_ids = [item for sublist in run_ids for item in sublist]
            data["run_id"] = list(set(run_ids))
        return data

    def get_pairing_from_argos_inputs(self, argos_inputs):
        sample_pairing = ""
        for i, job in enumerate(argos_inputs):
            tumor_sample_name = job["tumor"]["ID"]
            normal_sample_name = job["normal"]["ID"]
            sample_pairing += "\t".join([normal_sample_name, tumor_sample_name]) + "\n"
        return sample_pairing

    def get_mapping_from_argos_inputs(self, argos_inputs):
        sample_mapping = ""
        check_for_duplicates = list()
        filepaths = list()
        for i, job in enumerate(argos_inputs):
            tumor_sample_name = job["tumor"]["ID"]
            for p in job["tumor"]["R1"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([tumor_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["tumor"]["R2"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([tumor_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["tumor"]["zR1"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([tumor_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["tumor"]["zR2"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([tumor_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["tumor"]["bam"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([tumor_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)

            normal_sample_name = job["normal"]["ID"]
            for p in job["normal"]["R1"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([normal_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["normal"]["R2"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([normal_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["normal"]["zR1"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([normal_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["normal"]["zR2"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([normal_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
            for p in job["normal"]["bam"]:
                filepath = FileProcessor.parse_path_from_uri(p["location"])
                file_str = "\t".join([normal_sample_name, filepath]) + "\n"
                if file_str not in check_for_duplicates:
                    check_for_duplicates.append(file_str)
                    sample_mapping += file_str
                if filepath not in filepaths:
                    filepaths.append(filepath)
        return sample_mapping, filepaths

    def get_samples_from_patient_id(self, patient_id):
        """
        Retrieves samples from the database based on the patient_id

        Only retrieve patients from LIMS file group
        """
        patient_files = FileRepository.filter(
            file_group=os.environ.get("FILE_GROUP"),
            metadata=[f"cmoPatientId:{patient_id}"],
        )
        data = list()
        for current_file in patient_files:
            sample = dict()
            sample["id"] = current_file.file_id
            sample["path"] = current_file.path
            sample["file_name"] = current_file.file_name
            sample["metadata"] = current_file.metadata
            data.append(sample)

        samples = list()
        # group by igoId
        igo_id_group = dict()
        for sample in data:
            igo_id = sample["metadata"]["primaryId"]
            if igo_id not in igo_id_group:
                igo_id_group[igo_id] = list()
            igo_id_group[igo_id].append(sample)

        for igo_id in igo_id_group:
            samples.append(self.build_sample(igo_id_group[igo_id]))
        samples, bad_samples = self.remove_with_caveats(samples)
        number_of_bad_samples = len(bad_samples)
        if number_of_bad_samples > 0:
            LOGGER.warning(
                "Some samples for patient query %s have invalid %i values",
                patient_id,
                number_of_bad_samples,
            )
        return samples

    def compile_pairs(self, samples):
        """
        Creates pairs of tumors and normals from a list of samples
        """
        tumors = get_by_tumor_type(samples, "Tumor")
        normals = get_by_tumor_type(samples, "Normal")

        # pairing
        pairs = dict()
        pairs["tumor"] = list()
        pairs["normal"] = list()

        num_tumors = len(tumors)
        if num_tumors == 0:
            LOGGER.error("No tumor samples found; pairing will not be performed.")
            LOGGER.error("Returning an empty list of pairs.")

        for tumor in tumors:
            LOGGER.info("Pairing tumor sample %s", tumor["sample_id"])
            patient_id = tumor["patient_id"]
            if patient_id:
                run_mode = get_run_mode(tumor["run_mode"])
                bait_set = tumor["bait_set"]
                run_ids = tumor["run_id"]
                preservation_types = tumor["preservation_type"]
                normal = get_viable_normal(normals, patient_id, run_mode, bait_set)
                if normal:
                    LOGGER.info(
                        "Pairing %s (%s) with %s (%s)",
                        tumor["sample_id"],
                        tumor["SM"],
                        normal["sample_id"],
                        normal["SM"],
                    )
                    pairs["tumor"].append(tumor)
                    pairs["normal"].append(normal)
                else:
                    LOGGER.info(
                        "Missing normal for sample %s (%s); querying patient %s",
                        tumor["sample_id"],
                        tumor["SM"],
                        patient_id,
                    )
                    patient_samples = self.get_samples_from_patient_id(patient_id)
                    new_normals = get_by_tumor_type(patient_samples, "Normal")
                    new_normal = get_viable_normal(
                        new_normals, patient_id, run_mode, bait_set
                    )
                    if new_normal:
                        LOGGER.info(
                            "Pairing %s (%s) with %s (%s)",
                            tumor["sample_id"],
                            tumor["SM"],
                            new_normal["sample_id"],
                            new_normal["SM"],
                        )
                        pairs["tumor"].append(tumor)
                        pairs["normal"].append(new_normal)
                    else:
                        LOGGER.error(
                            "No normal found for %s (%s), patient %s",
                            tumor["sample_id"],
                            tumor["SM"],
                            patient_id,
                        )
        return pairs


class Fastqs:
    """
    Fastqs class to hold pairs of fastqs

    Does the pairing from a list of files

    The paired bool is True if all of the R1s in file list find a matching R2
    """

    def __init__(self, sample_name, file_list):
        self.sample_name = sample_name
        self.fastqs = dict()
        self.r1 = list()
        self.r2 = list()
        self.r2_bids = list()
        self.paired = True
        self._set_R(file_list)
        self.r1_bids = self._set_bids(self.r1)
        self.r2_bids = self._set_bids(self.r2)
        self.pu = self._set_pu()
        self.rg_id = self._set_rg_id()

    def _set_bids(self, r):
        r_bids = list()
        for f in r:
            r_file = get_file(f)
            r_bids.append(r_file.file_id)
        return r_bids

    def _set_pu(self):
        """
        Creating a list of PU values; used by argos pipeline as scatter input

        Only iterating across r1s since r1 and r2 should have the same metadata
        """
        pu = list()
        for f in self.r1:
            metadata = get_file(f).metadata
            flowcell_id = "MT_FCID"
            if "poolednormal" in self.sample_name.lower():
                flowcell_id = "PN_FCID"
                r = get_r_orientation(f)
                barcode_index = spoof_barcode(os.path.basename(f), r)
            else:
                fid = metadata["flowCellId"]
                if fid:
                    flowcell_id = fid
                barcode_index = metadata["barcodeIndex"]
            platform_unit = flowcell_id
            if barcode_index:
                platform_unit = "_".join([flowcell_id, barcode_index])
            pu.append(platform_unit)
        return pu

    def _set_rg_id(self):
        """
        Creating a list of RG_ID values; used by argos pipeline as scatter input

        Only iterating across r1s since r1 and r2 should have the same metadata
        """
        rg_ids = list()
        for i, f in enumerate(self.r1):
            metadata = get_file(f).metadata
            sample_name = self.sample_name
            pu = self.pu[i]
            rg_id = "_".join([sample_name, pu])
            rg_ids.append(rg_id)
        return rg_ids

    def _set_R(self, file_list):
        """
        From the file list, retrieve R1 and R2 fastq files

        Sets PU and bids, as well

        Uses _get_fastq_from_list() to find R2 pair.
        """
        r1s = list()
        r2s = list()
        for i in file_list:
            r = get_r_orientation(i.path)
            if r == "R1":
                r1s.append(i)
            if r == "R2":
                r2s.append(i)
        for f in r1s:
            self.r1.append(f.path)
            fastq1 = f.path
            expected_r2 = "R2".join(fastq1.rsplit("R1", 1))
            fastq2 = self._get_fastq_from_list(expected_r2, r2s)
            if fastq2:
                self.r2.append(fastq2.path)
            else:
                print("No fastq R2 found for %s" % f.path)
                self.paired = False

    def __str__(self):
        s = "R1:\n"
        for i in self.r1:
            s += i.path + "\n"
        s += "\nR2:\n"
        for i in self.r2:
            s += i.path + "\n"
        return s

    def _get_fastq_from_list(self, fastq_path, fastq_files):
        """
        Given fastq_path, find it in the list of fastq_files and return
        that File object
        """
        for f in fastq_files:
            fpath = f.path
            if fastq_path == fpath:
                return f


def get_file(fpath):
    data = FileRepository.filter(path=fpath)
    if data:
        return next(data)
    return None


def format_sample_name(sample_name, specimen_type, ignore_sample_formatting=False):
    """
    Formats a given sample_name to legacy ROSLIN naming conventions, provided that
    it is in valid CMO Sample Name format (see sample_pattern regex value, below)

    Current format is to prepend sample name with "s_" and convert all hyphens to
    underscores

    If it does not meet sample_pattern requirements OR is not a specimen_type=="CellLine",
    return 'sampleMalFormed'

    ignore_sample_formatting is applied if we want to return a sample name regardless of
    formatting
    """
    sample_pattern = re.compile(r"^[^0-9].*$")

    if not ignore_sample_formatting:
        try:
            if "s_" in sample_name[:2]:
                return sample_name
            elif (
                bool(sample_pattern.match(sample_name))
                or "cellline" in specimen_type.lower()
            ):  # cmoSampleName is formatted properly
                sample_name = "s_" + sample_name.replace("-", "_")
                return sample_name
            return sample_name
        except TypeError:
            LOGGER.error(
                "sampleNameError: sampleName is Nonetype; returning 'sampleNameMalformed'."
            )
            return "sampleNameMalformed"
    else:
        return sample_name


def get_run_mode(run_mode):
    if isinstance(run_mode, str):
        if "hiseq" in run_mode.lower():
            return "hiseq"
        if "novaseq" in run_mode.lower():
            return "novaseq"
        return run_mode


def get_r_orientation(fastq_filename):
    """
    Retrieve R orientation of fastq filename
    """
    reversed_filename = "".join(reversed(fastq_filename))
    r1_idx = reversed_filename.find("1R")
    r2_idx = reversed_filename.find("2R")
    if r1_idx == -1 and r2_idx == -1:
        return "ERROR"
    elif r1_idx > 0 and r2_idx == -1:
        return "R1"
    elif r2_idx > 0 and r1_idx == -1:
        return "R2"
    elif r1_idx > 0 and r2_idx > 0:
        if r1_idx < r2_idx:
            return "R1"
        return "R2"
    return "ERROR"


def get_by_tumor_type(data, tumor_type):
    """
    Retrieves a set of samples that contain a value tumor_type

    tumor_tupe is typically Normal or Tumor
    """
    samples = list()
    for sample in data:
        if tumor_type.lower() in sample["tumor_type"].lower():
            samples.append(sample)
    return samples


def get_viable_normal(normals, patient_id, run_mode, bait_set):
    """
    From a set of normals, return the ones that have matching patient_id, bait_set,
    run_mode, and the most recent

    Does not check for Pooled Normals; that's done separately
    """
    viable_normal = dict()
    for normal in normals:
        normal_run_mode = get_run_mode(normal["run_mode"])
        if (
            normal["patient_id"] == patient_id
            and normal["bait_set"] == bait_set
            and normal_run_mode == run_mode
        ):
            if viable_normal:
                try:
                    viable_normal = compare_dates(normal, viable_normal, "%y-%m-%d")
                except ValueError:
                    LOGGER.debug("Trying different date parser")
                    viable_normal = compare_dates(normal, viable_normal, "%Y-%m-%d")
            else:
                viable_normal = normal
    return viable_normal


def compare_dates(normal, viable_normal, date_string):
    """
    Compares dates between two normals; returns the most recent
    """
    for run_date in normal["run_date"]:
        normal_date = dt.strptime(run_date, date_string)
        for vrun_date in viable_normal["run_date"]:
            vnormal_date = dt.strptime(vrun_date, date_string)
            if vnormal_date < normal_date:
                viable_normal = normal
    return viable_normal


def spoof_barcode(sample_file_name, r_orientation):
    """
    Spoof barcode by removing 'R1' or 'R2' from the filename; paired fastqs
    are assumed to have only these two values as different

    We are also assuming there are no periods in the file names other than extensions
    """
    reversed_str = "".join(reversed(sample_file_name))
    if r_orientation == "R1":
        reversed_str = reversed_str.replace("1R", "")
    else:
        reversed_str = reversed_str.replace("2R", "")
    reversed_str = "".join(reversed(reversed_str))
    spoofed_barcode = reversed_str.split(os.extsep)[0]
    return spoofed_barcode


def normalize_igo_text_field(igo_text):
    # Flatten text data from the Genomics Core
    # to allow robust exact text matching.
    #
    # Allow variance in case and ignore non
    # alphanumeric characters (FYI).
    # Convert to lowercase
    s = igo_text.lower()
    # Remove special characters and extra spaces
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def format_sample(data, specimen_type):
    sample = dict()
    sample["ID"] = data["SM"]  # TODO: change someday
    sample["CN"] = data["CN"]
    sample["LB"] = data["LB"]
    sample["PL"] = data["PL"]
    sample["PU"] = data["PU"]
    sample["R1"] = list()
    sample["R2"] = list()
    sample["zR1"] = list()  # TODO: Add for Xenografts
    sample["zR2"] = list()  # TODO: Add for Xenografts
    sample["bam"] = list()
    sample["zBam"] = list()
    sample["RG_ID"] = data["ID"]
    sample["adapter"] = (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG"
    )
    sample["adapter2"] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    sample["bwa_output"] = sample["ID"] + ".bam"
    sample["request_id"] = data["request_id"]
    sample["specimen_type"] = data["specimen_type"]

    if specimen_type in PDX_SPECIMEN_TYPES:
        r1 = "zR1"
        r2 = "zR2"
        bam = "zBam"
    elif specimen_type in NON_PDX_SPECIMEN_TYPES:
        r1 = "R1"
        r2 = "R2"
        bam = "bam"
    else:
        raise Exception(f"Invalid Specimen Type: {specimen_type}")

    for i in data["R1"]:
        if i:
            sample[r1].append({"class": "File", "location": "juno://" + i})
    for i in data["R2"]:
        if i:
            sample[r2].append({"class": "File", "location": "juno://" + i})
    for i in data["bam"]:
        if i:
            sample[bam].append({"class": "File", "location": "juno://" + i})
    return sample


def convert_references(project_id, assay, pi, pi_email):
    genomic_resources = load_references()
    request_files = genomic_resources["request_files"]
    intervals = get_baits_and_targets(assay, genomic_resources)
    curated_bams = get_curated_bams(assay, request_files)
    covariates = [
        "CycleCovariate",
        "ContextCovariate",
        "ReadGroupCovariate",
        "QualityScoreCovariate",
    ]
    rf = ["BadCigar"]
    genome = "GRCh37"
    delly_type = ["DUP", "DEL", "INV", "INS", "BND"]
    facets_cval = get_facets_cval(assay)
    facets_pcval = get_facets_pcval(assay)
    complex_nn = get_complex_nn(assay)
    complex_tn = get_complex_tn(assay)
    temp_dir = "/scratch"
    if "TMPDIR" in os.environ:
        if os.environ["TMPDIR"]:
            temp_dir = os.environ["TMPDIR"]

    files = {
        "refseq": {"class": "File", "location": str(request_files["refseq"])},
        "vep_data": str(request_files["vep_data"]),
        "hotspot_list": str(request_files["hotspot_list"]),
        "hotspot_list_maf": {
            "class": "File",
            "location": str(request_files["hotspot_list_maf"]),
        },
        "delly_exclude": {
            "class": "File",
            "location": str(genomic_resources["genomes"][genome]["delly"]),
        },
        "hotspot_vcf": str(request_files["hotspot_vcf"]),
        "facets_snps": {
            "class": "File",
            "location": str(genomic_resources["genomes"][genome]["facets_snps"]),
        },
        "custom_enst": str(request_files["custom_enst"]),
        "vep_path": str(request_files["vep_path"]),
        "conpair_markers": str(request_files["conpair_markers"]),
        "conpair_markers_bed": str(request_files["conpair_markers_bed"]),
    }

    files.update(intervals)

    out_dict = {
        "curated_bams": curated_bams,
        "hapmap": {"class": "File", "location": str(request_files["hapmap"])},
        "dbsnp": {"class": "File", "location": str(request_files["dbsnp"])},
        "indels_1000g": {
            "class": "File",
            "location": str(request_files["indels_1000g"]),
        },
        "snps_1000g": {"class": "File", "location": str(request_files["snps_1000g"])},
        "cosmic": {"class": "File", "location": str(request_files["cosmic"])},
        "exac_filter": {"class": "File", "location": str(request_files["exac_filter"])},
        "ref_fasta": {"class": "File", "location": str(request_files["ref_fasta"])},
        "mouse_fasta": {"class": "File", "location": str(request_files["mouse_fasta"])},
        "db_files": files,
    }
    # emit_original_quals boolean could be problematic; test
    params = {
        "abra_scratch": temp_dir,
        "abra_ram_min": 84000,
        "genome": genome,
        "intervals": genomic_resources["genomes"][genome]["intervals"],
        "mutect_dcov": 50000,
        "mutect_rf": rf,
        "num_cpu_threads_per_data_thread": 6,
        "covariates": covariates,
        "emit_original_quals": True,
        "num_threads": 10,
        "assay": assay,
        "tmp_dir": temp_dir,
        "project_prefix": project_id,
        "opt_dup_pix_dist": "2500",
        "delly_type": delly_type,
        "facets_cval": facets_cval,
        "facets_pcval": facets_pcval,
        "complex_nn": complex_nn,
        "complex_tn": complex_tn,
        "scripts_bin": "/usr/bin",
        "gatk_jar_path": "/usr/bin/gatk.jar",
        "pi": pi,
        "pi_email": pi_email,
    }
    out_dict.update({"runparams": params})
    return out_dict


def load_references():
    WORKDIR = os.path.dirname(os.path.abspath(__file__))
    d = json.load(
        open(os.path.join(WORKDIR, "reference_jsons/genomic_resources.json"), "rb")
    )
    return d


def get_baits_and_targets(assay, genomic_resources):
    # probably need similar rules for whatever "Exome" string is in rquest
    targets = genomic_resources["targets"]

    target_assay = assay

    if assay.find("HemePACT_v4") > -1:
        target_assay = "HemePACT_v4_BAITS"

    if assay.find("IMPACT-Heme_v2") > -1:
        target_assay = "IMPACT-Heme_v2_BAITS"

    if assay.find("IMPACT505") > -1:
        target_assay = "IMPACT505_b37"
    if assay.find("IMPACT410") > -1:
        target_assay = "IMPACT410_b37"
    if assay.find("IMPACT468") > -1:
        target_assay = "IMPACT468_b37"
    if assay.find("IMPACT341") > -1:
        target_assay = "IMPACT341_b37"
    if assay.find("IDT_Exome_v1_FP") > -1:
        target_assay = "IDT_Exome_v1_FP_b37"
    if assay.find("IMPACT468+08390") > -1:
        target_assay = "IMPACT468_08390"
    if assay.find("IMPACT468+Poirier_RB1_intron_V2") > -1:
        target_assay = "IMPACT468_08050"

    if target_assay in targets:
        return {
            "bait_intervals": {
                "class": "File",
                "location": str(targets[target_assay]["baits_list"]),
            },
            "target_intervals": {
                "class": "File",
                "location": str(targets[target_assay]["targets_list"]),
            },
            "fp_intervals": {
                "class": "File",
                "location": str(targets[target_assay]["FP_intervals"]),
            },
            "fp_genotypes": {
                "class": "File",
                "location": str(targets[target_assay]["FP_genotypes"]),
            },
        }
    else:
        LOGGER.error(
            "ERROR: Targets for Assay not found in genomic_resources.json: %s", assay
        )


def get_facets_cval(assay):
    if assay.find("IMPACT") > -1 or assay.find("HemePACT") > -1:
        return 50
    return 100


def get_facets_pcval(assay):
    if assay.find("IMPACT") > -1 or assay.find("HemePACT") > -1:
        return 100
    return 500


def get_complex_nn(assay):
    if assay.find("IMPACT") > -1 or assay.find("HemePACT") > -1:
        return 0.2
    return 0.1


def get_curated_bams(assay, request_files):
    # Default to AgilentExon_51MB_b37_v3 BAMs for all assays except those specified below
    json_curated_bams = request_files["curated_bams"]["AgilentExon_51MB_b37_v3"]

    if assay.find("IMPACT-Heme") > -1:
        json_curated_bams = request_files["curated_bams"]["IMPACT-Heme_v2_BAITS"]
    # Default to IMPACT468_b37 BAMs for all IMPACT/HemePACT assays
    elif assay.find("IMPACT") > -1 or assay.find("HemePACT") > -1:
        json_curated_bams = request_files["curated_bams"]["IMPACT468_b37"]
    elif assay.find("IDT_Exome_v1_FP") > -1:
        json_curated_bams = request_files["curated_bams"]["IDT_Exome_v1_FP_b37"]
    array = []
    for bam in json_curated_bams:
        array.append({"class": "File", "location": str(bam)})
    return array


def get_complex_tn(assay):
    if assay.find("IMPACT") > -1 or assay.find("HemePACT") > -1:
        return 0.5
    return 0.2
