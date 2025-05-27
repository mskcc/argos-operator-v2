import sys
import pprint

sys.path.append("..")
import argos_operator_v2.argos_operator_v2
from argos_operator_v2.argos_operator_v2 import ArgosOperatorV2

pipeline = {
        "pipeline_id": "7af9f6c2-3820-11f0-b5d5-ac1f6bb4ad16",
        "pipeline_name": "incl_bedfile",
        "pipeline_link": "https://github.com/mskcc/argos-cwl",
        "pipeline_version": "incl_bedfile",
        "pipeline_entrypoint": "workflows/pair-workflow-sv.cwl",
        "pipeline_format": 0
    }

op = ArgosOperatorV2(request_id="08944_B", pipeline=pipeline)

run_ids = ["FAUCI2_BLAH"]
bait_set = "IMPACT505_BAITS"
preservation_types = ["FROZEN"]

response = op.get_pooled_normal(run_ids=run_ids, bait_set=bait_set, preservation_types=preservation_types)

pprint.pprint(response)
