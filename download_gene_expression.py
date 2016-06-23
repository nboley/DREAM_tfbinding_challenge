import os, sys

from pyTFbindtools.ENCODE_ChIPseq_tools import (
    find_ENCODE_RNAseq_experiment_ids,
    find_ENCODE_DCC_experiment_files)
from pyTFbindtools.DB import conn

def load_samples():
    curr = conn.cursor()
    query = "select * from dream_challenge.challenge_sample_types;"
    curr.execute(query)
    return set(x[0] for x in curr.fetchall())

def load_eid_mapping():
    curr = conn.cursor()
    query = "select sample_type, eid from eid_mapping;"
    curr.execute(query)
    return dict(curr.fetchall())

def main():
    sample_types = load_samples()
    eid_mapping = load_eid_mapping()
    for exp_id in find_ENCODE_RNAseq_experiment_ids('hg19'):
        files = list(find_ENCODE_DCC_experiment_files(exp_id))
        if files[0].sample_type not in sample_types: continue
        any_hg19 = any(f.assembly == 'hg19' for f in files)
        for f in files:
            if f.output_type != 'gene quantifications': continue
            if f.file_format != 'tsv': continue
            if any_hg19 and f.assembly != 'hg19':
                continue
            try: eid = eid_mapping[f.sample_type]
            except KeyError: eid = "E999"
            ofname = "gene_expression.{eid}.{sample_type}.{0.bsid}.{0.exp_id}.tsv".format(
                f,
                eid=eid,
                sample_type=f.sample_type.replace(" ", "_"))
            print f.sample_type, ofname
            os.system("wget https://encodeproject.org/{} -O {}".format(f.file_loc, ofname)) 
        print exp_id
    return

if __name__ == '__main__':
    main()
