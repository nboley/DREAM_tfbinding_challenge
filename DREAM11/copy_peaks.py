import os, sys

from subprocess import check_output
from collections import namedtuple, defaultdict

import numpy as np

import pybedtools

from grit.lib.multiprocessing_utils import run_in_parallel

from build_region_labels import build_regions_labels_from_beds

METADATA_TSV = "/mnt/data/TF_binding/DREAM_challenge/DREAM_TFs_FINAL_SAMPLE_SHEET.v3.tsv"

MetaDataRecord = namedtuple('MetaDataRecord', [
    "IDX", "SAMPLE_NAME", "TF", "CELL_TYPE", "ENCODE_ID", "HIDDEN_TEST_SET",
    "TODO", "NOT_RELEASED", "LOWEST_QUALITY_REPLICATE", "MISSING_RNASEQ"])

samples_to_skip = []
#"CHIPseq.E2F6.HeLa-S3.EXPID_ENCSR000EVK"

TRAIN_REGIONS_BED = "/mnt/data/TF_binding/DREAM_challenge/public_data/annotations/train_regions.blacklisted.bed.gz"
#TRAIN_REGIONS_BED = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/labels/tmp.bed.gz"
CHIPSEQ_IDR_PEAKS_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/idr/"
CHIPSEQ_RELAXED_PEAKS_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/naive_overlap/"
CHIPSEQ_FOLD_CHANGE_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/fold_change_signal/"
CHIPSEQ_LABELS = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/labels/"
DNASE_IDR_PEAKS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/idr/"


def find_num_peaks(fname):
    return int(check_output("zcat {} | wc -l".format(fname) , shell=True))

def sample_and_factor_from_samplename(sample_name):
    return tuple(sample_name.split(".")[1:3])

def load_metadata(fname=METADATA_TSV):
    all_data = []
    with open(fname) as fp:
        for i, line in enumerate(fp):
            # skip the header
            if i == 0: continue
            data = line.strip().split("\t")
            data[0] = int(data[0])
            data[-5:] = [bool(int(x)) for x in data[-5:]]
            all_data.append(MetaDataRecord(*data))
    return all_data

def find_idr_peaks(sample_names, base_dir, suffix):
    for sample_name in sample_names:
        if sample_name in samples_to_skip: continue
        try:
            for fname in os.listdir(base_dir.format(sample_name)):
                if not fname.endswith(suffix): continue
                yield sample_name, os.path.join(base_dir.format(sample_name), fname)
        except OSError, inst:
            print inst
            continue
    return

def copy_chipseq_data(
        sample_dirs,
        idr_peaks_output_dir,
        relaxed_peaks_output_dir,
        fold_change_output_dir,
        regions_bed_fname):
    ## first group the IDR optimal peaks, and when there are alternates choose
    ## the experiment that has the largest number of peaks
    optimal_peaks = defaultdict(list)
    for sample_name, pk_fname in find_idr_peaks(
           sample_dirs,
           "/mnt/data/TF_binding/DREAM_challenge/all_data/chipseq/output/DREAM_challenge/{}/out/peak/idr/optimal_set/",
           ".IDR0.05.filt.narrowPeak.gz"):
        optimal_peaks[sample_and_factor_from_samplename(sample_name)].append(
            (sample_name, pk_fname))
    # the new sample_dirs has a single entry for each TF,sample_type combo
    sample_dirs = []
    for (factor, sample), sample_and_fnames in optimal_peaks.iteritems():
        if len(sample_and_fnames) == 1:
            sample_dirs.append(sample_and_fnames[0][0])
        elif len(sample_and_fnames) > 1:
            assert False, "There shouldn't be any samples with more than 1 replicate"
            #print sample_and_fnames[0][0], len(sample_and_fnames)
            best_sample, most_num_lines = None, 0
            for sample_name, fname in sample_and_fnames:
                num_lines = find_num_peaks(fname)
                if num_lines > most_num_lines:
                    best_sample = sample_name
                    most_num_lines = num_lines
            print best_sample, most_num_lines
        else:
            assert False, "It shouldn't be possible to have zero fnames"
    ## now that we've refined the sample list, copy all of the peaks and wiggles
    ## into the correct directory
    # copy the IDR peaks
    cmds = []
    for sample_name, pk_fname in find_idr_peaks(
           sample_dirs,
           "/mnt/data/TF_binding/DREAM_challenge/all_data/chipseq/output/DREAM_challenge/{}/out/peak/idr/optimal_set/",
           ".IDR0.05.filt.narrowPeak.gz"):
        sample, factor = sample_and_factor_from_samplename(sample_name)
        ofname = "ChIPseq.{factor}.{sample}.conservative.train.narrowPeak.gz".format(
            factor=factor, sample=sample)
        cmd = "bedtools intersect -wa -u -a {pk_fname} -b {regions_fname} | pigz > {ofname} ".format(
            pk_fname=pk_fname,
            regions_fname=regions_bed_fname,
            ofname=os.path.join(idr_peaks_output_dir, ofname))
        cmds.append([cmd,])
    # copy the relaxed peaks
    for sample_name, pk_fname in find_idr_peaks(
            sample_dirs,
            "/mnt/data/TF_binding/DREAM_challenge/all_data/chipseq/output/DREAM_challenge/{}/out/peak/idr/pooled_pseudo_reps/",
            "unthresholded-peaks.txt.gz"):
        sample, factor = sample_and_factor_from_samplename(sample_name)
        ofname = "ChIPseq.{factor}.{sample}.relaxed.narrowPeak.gz".format(
            factor=factor, sample=sample)
        cmd = "bedtools intersect -wa -u -a {pk_fname} -b {regions_fname} | pigz > {ofname} ".format(
            pk_fname=pk_fname,
            regions_fname=regions_bed_fname,
            ofname=os.path.join(relaxed_peaks_output_dir, ofname))
        cmds.append([cmd,])

    # run all of the peaks intersection cmds
    run_in_parallel(16, os.system, cmds)

    # copy the fold change wiggles
    fold_change_bigwigs = list(find_idr_peaks(
        sample_dirs,
        "/mnt/data/TF_binding/DREAM_challenge/all_data/chipseq/output/DREAM_challenge/{}/out/signal/macs2/pooled_rep/",
        ".fc.signal.bw"))
    cmds = []
    for i, (sample_name, pk_fname) in enumerate(fold_change_bigwigs):
        #print i, len(fold_change_bigwigs)
        try:
            sample, factor = sample_and_factor_from_samplename(sample_name)
            ofname = "ChIPseq.{factor}.{sample}.fc.signal.train.bw".format(
                factor=factor, sample=sample)
            try:
                open(ofname)
            except IOError:
                cmd = "bwtool remove mask -inverse {} {} {}".format(
                    regions_bed_fname,
                    pk_fname,
                    os.path.join(fold_change_output_dir, ofname))
                cmds.append([cmd,])
            else:
                pass
            ### the old copy command
            #os.system("cp -u {} {}".format(
            #    pk_fname, os.path.join(fold_change_output_dir, ofname)))
        except FileNotFoundError:
            print "Can't run:", cmd
    run_in_parallel(16, os.system, cmds)

def copy_public_chipseq_data():
    metadata = load_metadata()
    train_sample_dirs = set(
        x.SAMPLE_NAME for x in metadata
        if (x.SAMPLE_NAME != ''
            and not x.TODO
            and not x.HIDDEN_TEST_SET
            and not x.LOWEST_QUALITY_REPLICATE
            and not x.MISSING_RNASEQ)
    )
    copy_chipseq_data(
        sample_dirs=train_sample_dirs,
        idr_peaks_output_dir=CHIPSEQ_IDR_PEAKS_DIR,
        relaxed_peaks_output_dir=CHIPSEQ_RELAXED_PEAKS_DIR,
        fold_change_output_dir=CHIPSEQ_FOLD_CHANGE_DIR,
        regions_bed_fname="/mnt/data/TF_binding/DREAM_challenge/public_data/annotations/train_regions.blacklisted.merged.bed"
    )
    pass

def copy_private_chipseq_data():
    metadata = load_metadata()
    test_sample_dirs = set(
        x.SAMPLE_NAME for x in metadata
        if (x.SAMPLE_NAME != ''
            and not x.TODO
            and x.HIDDEN_TEST_SET
            and not x.LOWEST_QUALITY_REPLICATE
            and not x.MISSING_RNASEQ)
    )
    copy_chipseq_data(
        sample_dirs=test_sample_dirs,
        idr_peaks_output_dir="/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/idr/",
        relaxed_peaks_output_dir="/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/naive_overlap/",
        fold_change_output_dir="/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/fold_change_signal/"
    )
    return



def copy_DNASE_files():
    # find which experiment IDS to copy
    grpd_peak_files = defaultdict(list)
    for fname in os.listdir(DNASE_IDR_PEAKS_PATH):
        sample_type = fname.split(".")[1].split("_")[0]
        grpd_peak_files[sample_type].append(os.path.join(DNASE_IDR_PEAKS_PATH, fname))

    sample_types = []
    for sample_type, fnames in grpd_peak_files.iteritems():
        if len(fnames) == 1:
            sample_types.append(sample_type)
        elif len(fnames) > 1:
            assert False
            print sample_type, len(fnames)
            for fname in fnames:
                bsid = 1 if "b1.36bp" in fname else 2
                print find_num_peaks(fname), fname

    # copy the IDR peaks
    DNASE_IDR_PEAKS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/idr/"
    for fname in os.listdir(DNASE_IDR_PEAKS_PATH):
        sample_type = fname.split(".")[1].split("_")[0]
        ofname = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/idr/DNASE.{sample_type}.conservative.narrowPeak.gz"
        cmd = "cp -u {} {}".format(os.path.join(DNASE_IDR_PEAKS_PATH, fname),
                                ofname.format(sample_type=sample_type))
        print cmd
        os.system(cmd)
        
    # copy the relaxed peaks
    DNASE_RELAXED_PEAKS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/naive_overlap/"
    for fname in os.listdir(DNASE_RELAXED_PEAKS_PATH):
        sample_type = fname.split(".")[1].split("_")[0]
        ofname = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/naive_overlap/DNASE.{sample_type}.relaxed.narrowPeak.gz"
        cmd = "cp -u {} {}".format(os.path.join(DNASE_RELAXED_PEAKS_PATH, fname),
                                ofname.format(sample_type=sample_type))
        print cmd
        os.system(cmd)

    # copy the bam files
    bam_counter = defaultdict(int)
    DNASE_BAMS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/bams/"
    for fname in os.listdir(DNASE_BAMS_PATH):
        sample_type = fname.split(".")[1].split("_")[0]
        if '.b1.' in fname: biorep = 1
        elif '.b2.' in fname: biorep = 2
        else: assert False, fname
        bam_counter[(sample_type, biorep)] += 1
        ofname = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/bams/DNASE.{sample_type}.biorep{biorep}.techrep{rep_id}.bam"
        cmd = "cp -u {} {}".format(
            os.path.join(DNASE_BAMS_PATH, fname),
            ofname.format(sample_type=sample_type,
                          biorep=biorep,
                          rep_id=bam_counter[(sample_type, biorep)])
        )
        print cmd
        os.system(cmd)
    
    pass

def build_train_region_wiggles():
    pass

def build_labels_for_sample_and_factor(
        regions_fname,
        idr_peaks_fname, relaxed_peaks_fname,
        output_directory,
        overwrite=False):
    sample, factor = idr_peaks_fname.split(".")[1:3]
    assert [sample, factor] == relaxed_peaks_fname.split(".")[1:3]
    ofname = os.path.join(
        output_directory, "%s.%s.train.labels" % (sample, factor))
    if not overwrite:
        try:
            open(ofname + ".npy")
        except IOError:
            pass
        else:
            print "Skipping %s-%s: labels matrix file already exists" % (
                sample, factor)
            return
    print "Processing ", sample, factor
    regions_bed = pybedtools.BedTool(regions_fname) #.sort() - already sorted
    relaxed_peaks_bed = pybedtools.BedTool(relaxed_peaks_fname).sort().merge()
    idr_peaks_bed = pybedtools.BedTool(idr_peaks_fname).sort().merge()

    labels = build_regions_labels_from_beds(
        regions_bed, idr_peaks_bed, relaxed_peaks_bed)
    print "Saving to %s" % ofname
    np.save(ofname, labels)
    return

def build_labels(overwrite_existing=False):
    matched_peaks = defaultdict(lambda: {'idr': None, 'relaxed': None})
    for fname in os.listdir(CHIPSEQ_IDR_PEAKS_DIR):
        sample, factor = fname.split(".")[1:3]
        assert matched_peaks[(sample, factor)]['idr'] is None
        matched_peaks[(sample, factor)]['idr'] = os.path.join(
            CHIPSEQ_IDR_PEAKS_DIR, fname)

    for fname in os.listdir(CHIPSEQ_RELAXED_PEAKS_DIR):
        sample, factor = fname.split(".")[1:3]
        assert matched_peaks[(sample, factor)]['idr'] is not None
        assert matched_peaks[(sample, factor)]['relaxed'] is None
        matched_peaks[(sample, factor)]['relaxed'] = os.path.join(
            CHIPSEQ_RELAXED_PEAKS_DIR, fname)

    all_args = []
    output_dir = \
        "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/labels/"
    for (sample, factor), fnames in matched_peaks.iteritems():
        all_args.append(
            (TRAIN_REGIONS_BED, fnames['idr'], fnames['relaxed'], output_dir))
    run_in_parallel(16, build_labels_for_sample_and_factor, all_args)
    return
    
def main():
    #copy_private_chipseq_data()
    copy_public_chipseq_data()
    #copy_DNASE_files()
    #build_labels()

if __name__ == '__main__':
    main()
