import os, sys

from subprocess import check_output
from collections import namedtuple, defaultdict

METADATA_TSV = "/mnt/data/TF_binding/DREAM_challenge/chipseq/peaks/DREAM_TFs_FINAL_SAMPLE_SHEET.v2.tsv"

MetaDataRecord = namedtuple('MetaDataRecord', [
    "IDX", "SAMPLE_NAME", "TF", "CELL_TYPE", "ENCODE_ID", "HIDDEN_TEST_SET", "TODO", "NOT_RELEASED", "LOWEST_QUALITY_REPLICATE"])

samples_to_skip = []
#"CHIPseq.E2F6.HeLa-S3.EXPID_ENCSR000EVK"

CHIPSEQ_IDR_PEAKS_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/idr/"
CHIPSEQ_RELAXED_PEAKS_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/naive_overlap/"
CHIPSEQ_FOLD_CHANGE_DIR = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/fold_change_signal/"

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
            data[-4:] = [bool(int(x)) for x in data[-4:]]
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

def copy_chipseq_data():
    metadata = load_metadata()
    train_sample_dirs = set(
        x.SAMPLE_NAME for x in metadata
        if (x.SAMPLE_NAME != ''
            and not x.TODO
            and not x.HIDDEN_TEST_SET
            and not x.LOWEST_QUALITY_REPLICATE)
    )
    
    ## first group the IDR optimal peaks, and when there are alternates choose
    ## the experiment that has the largest number of peaks
    optimal_peaks = defaultdict(list)
    for sample_name, pk_fname in find_idr_peaks(
           train_sample_dirs,
           "/mnt/data/TF_binding/DREAM_challenge/chipseq/output/DREAM_challenge/{}/out/peak/idr/optimal_set/",
           ".IDR0.05.filt.narrowPeak.gz"):
        optimal_peaks[sample_and_factor_from_samplename(sample_name)].append(
            (sample_name, pk_fname))
    # the new train_sample_dirs has a single entry for each TF,sample_type combo
    train_sample_dirs = []
    for (factor, sample), sample_and_fnames in optimal_peaks.iteritems():
        if len(sample_and_fnames) == 1:
            train_sample_dirs.append(sample_and_fnames[0][0])
        elif len(sample_and_fnames) > 1:
            assert False, "There shouldn't be any samples with more than 1 replicate"
            print sample_and_fnames[0][0], len(sample_and_fnames)
            best_sample, most_num_lines = None, 0
            for sample_name, fname in sample_and_fnames:
                num_lines = find_num_peaks(fname)
                if num_lines > most_num_lines:
                    best_sample = sample_name
                    most_num_lines = num_lines
            #print best_sample, most_num_lines
        else:
            assert False, "It shouldn't be possible to have zero fnames"

    ## now that we've refined the sample list, copy all of the peaks and wiggles
    ## into the correct directory
    # copy the IDR peaks
    for sample_name, pk_fname in find_idr_peaks(
           train_sample_dirs,
           "/mnt/data/TF_binding/DREAM_challenge/chipseq/output/DREAM_challenge/{}/out/peak/idr/optimal_set/",
           ".IDR0.05.filt.narrowPeak.gz"):
        sample, factor = sample_and_factor_from_samplename(sample_name)
        ofname = "ChIPseq.{factor}.{sample}.conservative.narrowPeak.gz".format(
            factor=factor, sample=sample)
        cmd = "cp {} {}".format(
            pk_fname, os.path.join(CHIPSEQ_IDR_PEAKS_DIR, ofname))
        os.system(cmd)
    # copy the relaxed peaks
    for sample_name, pk_fname in find_idr_peaks(
            train_sample_dirs,
            "/mnt/data/TF_binding/DREAM_challenge/chipseq/output/DREAM_challenge/{}/out/peak/idr/pooled_pseudo_reps/",
            "unthresholded-peaks.txt.gz"):
        sample, factor = sample_and_factor_from_samplename(sample_name)
        ofname = "ChIPseq.{factor}.{sample}.relaxed.narrowPeak.gz".format(
            factor=factor, sample=sample)
        os.system("cp {} {}".format(
            pk_fname, os.path.join(CHIPSEQ_RELAXED_PEAKS_DIR, ofname)))
    # copy the fold change wiggles
    fold_change_bigwigs = list(find_idr_peaks(
        train_sample_dirs,
        "/mnt/data/TF_binding/DREAM_challenge/chipseq/output/DREAM_challenge/{}/out/signal/macs2/pooled_rep/",
        ".fc.signal.bw"))
    for i, (sample_name, pk_fname) in enumerate(fold_change_bigwigs):
        print i, len(fold_change_bigwigs)
        try:
            sample, factor = sample_and_factor_from_samplename(sample_name)
            ofname = "ChIPseq.{factor}.{sample}.fc.signal.bw".format(
                factor=factor, sample=sample)
            os.system("cp {} {}".format(
                pk_fname, os.path.join(CHIPSEQ_FOLD_CHANGE_DIR, ofname)))
        except FileNotFoundError:
            print "Can't run:", cmd

def copy_DNASE_files():
    # find which experiment IDS to copy
    DNASE_IDR_PEAKS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/idr/"
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
        cmd = "cp {} {}".format(os.path.join(DNASE_IDR_PEAKS_PATH, fname),
                                ofname.format(sample_type=sample_type))
        print cmd
        os.system(cmd)
        
    # copy the relaxed peaks
    DNASE_RELAXED_PEAKS_PATH = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/naive_overlap/"
    for fname in os.listdir(DNASE_RELAXED_PEAKS_PATH):
        sample_type = fname.split(".")[1].split("_")[0]
        ofname = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/naive_overlap/DNASE.{sample_type}.relaxed.narrowPeak.gz"
        cmd = "cp {} {}".format(os.path.join(DNASE_RELAXED_PEAKS_PATH, fname),
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
        cmd = "cp {} {}".format(
            os.path.join(DNASE_BAMS_PATH, fname),
            ofname.format(sample_type=sample_type,
                          biorep=biorep,
                          rep_id=bam_counter[(sample_type, biorep)])
        )
        print cmd
        os.system(cmd)
    
    pass

def main():
    #copy_chipseq_data()
    copy_DNASE_files()

if __name__ == '__main__':
    main()
