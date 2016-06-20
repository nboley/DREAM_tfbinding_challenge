import os, sys
from collections import namedtuple

import pyTFbindtools.ENCODE_ChIPseq_tools
pyTFbindtools.ENCODE_ChIPseq_tools.BASE_URL = "https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/"

from pyTFbindtools.ENCODE_ChIPseq_tools import *

TFData = namedtuple('TFData', [
    'tf_name', 'sample_type', 'chipseq_fastqs', 'input_fastqs'])

def load_tfnames(fname="/mnt/data/TF_binding/DREAM_challenge/chipseq/factor_names"):
    factor_names = set()
    with open(fname) as fp:
        for line in fp:
            factor_names.add(line.strip())
    return factor_names

def build_anshuls_list(annotation):
    from pyTFbindtools.DB import conn
    factor_names = load_tfnames()
    unobserved_factornames = load_tfnames()
    cur = conn.cursor()
    query = """
    SELECT tf_name, 
           sample_type, 
           encode_experiment_id 
      FROM encode_chipseq_experiments_with_dnase NATURAL JOIN tfs;
    """
    cur.execute(query, [])
    res = cur.fetchall()
    for i, (tf_name,
            sample_type, 
            chipseq_exp_id) in enumerate(res):
        #print >> sys.stderr, i, len(res), chipseq_exp_id
        if tf_name not in factor_names: continue
        if tf_name in unobserved_factornames:
            unobserved_factornames.remove(tf_name)
        input_exp_ids = find_ENCODE_DCC_experiment_controls(chipseq_exp_id)
        chipseq_fastqs = find_ENCODE_DCC_fastqs(chipseq_exp_id)
        input_fastqs = [ find_ENCODE_DCC_fastqs(exp_id)
                         for exp_id in input_exp_ids ]
        data = TFData(tf_name, sample_type, chipseq_fastqs, input_fastqs)
        yield data

    print >> sys.stderr, unobserved_factornames

def load_tf_data(exp_id):
    input_exp_ids = find_ENCODE_DCC_experiment_controls(chipseq_exp_id)
    chipseq_fastqs = find_ENCODE_DCC_fastqs(chipseq_exp_id)
    input_fastqs = [ find_ENCODE_DCC_fastqs(exp_id)
                     for exp_id in input_exp_ids ]
    print chipseq_fastqs
    assert False
    data = TFData(tf_name, sample_type, chipseq_fastqs, input_fastqs)

    pass

def download_fastqs(data):
    assert len(data.chipseq_fastqs[1]) == len(data.chipseq_fastqs[2])
    for x in data.chipseq_fastqs[0]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.unpaired.fastq.gz" % (
                data.tf_name, data.sample_type.replace(" ", "_"))
        ).format(x)
        print "wget --quiet {} -O {}".format(BASE_URL + x.file_loc, ofname)
    for x in data.chipseq_fastqs[1]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.R1.fastq.gz" % (
                data.tf_name, data.sample_type.replace(" ", "_"))
        ).format(x)
        print "wget --quiet {} -O {}".format(BASE_URL + x.file_loc, ofname)
    for x in data.chipseq_fastqs[2]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.R2.fastq.gz" % (
                data.tf_name, data.sample_type.replace(" ", "_"))
        ).format(x)
        print "wget --quiet {} -O {}".format(BASE_URL + x.file_loc, ofname)
    return

def group_and_cat_controls():
    unpaired_controls = defaultdict(set)
    paired_controls = defaultdict(set)
    for i, data in enumerate(build_anshuls_list('hg19')):
        # deal with the unpaired files
        for control_exp in data.input_fastqs:
            # unpaired
            for control in control_exp[0]:
                unpaired_controls[control.sample_type].add(control)
            # paired 
            for control in zip(control_exp[1], control_exp[2]):
                assert control[0].sample_type == control[1].sample_type
                paired_controls[control[0].sample_type].add(
                    (control[0], control[1])
                )

    for sample_type, controls in unpaired_controls.iteritems():
        cat_files = []
        for i, x in enumerate(controls):
            ofname = "CONTROL.{0}.{1}.FILEID{2}.unpaired.fastq.gz".format(
                sample_type.replace(" ", "_").replace("/", "_"),
                str(i).zfill(3),
                x.id
            )
            print "wget   {} -O {}".format(BASE_URL + x.file_loc, ofname)
            cat_files.append(ofname)
        print "cat {} > CONTROL.{}.unpaired.fastq.gz".format(
            " ".join(cat_files), sample_type.replace(" ", "_").replace("/", "_") )
        print "rm {}".format(" ".join(cat_files))
    
    for sample_type, controls in paired_controls.iteritems():
        r1_cat_files = []
        r2_cat_files = []
        for i, (r1, r2) in enumerate(controls):
            ofname = "CONTROL.{0}.{1}.FILEID{2}.R1.fastq.gz".format(
                sample_type.replace(" ", "_").replace("/", "_"),
                str(i).zfill(3),
                r1.id
            )
            print "wget   {} -O {}".format(BASE_URL + r1.file_loc, ofname)
            r1_cat_files.append(ofname)
            ofname = "CONTROL.{0}.{1}.FILEID{2}.R2.fastq.gz".format(
                sample_type.replace(" ", "_").replace("/", "_"),
                str(i).zfill(3),
                r2.id
            )
            print "wget   {} -O {}".format(BASE_URL + r2.file_loc, ofname)
            r2_cat_files.append(ofname)
        print "cat {} > CONTROL.{}.R1.fastq.gz".format(
            " ".join(r1_cat_files), sample_type.replace(" ", "_").replace("/", "_") )
        print "rm {}".format(" ".join(r1_cat_files))
        print "cat {} > CONTROL.{}.R2.fastq.gz".format(
            " ".join(r2_cat_files), sample_type.replace(" ", "_").replace("/", "_") )
        print "rm {}".format(" ".join(r2_cat_files))

def download_public_data():
    anshuls_list = list(build_anshuls_list('hg19'))
    for i, data in enumerate(anshuls_list):
        print "echo ", i, len(anshuls_list)
        print >> sys.stderr, i, len(anshuls_list), data.tf_name, data.sample_type
        download_fastqs(data)

def main():
    with open(sys.argv[1]) as fp:
        pass
    load_tf_data(exp_id)
if __name__ == '__main__':
    main()
