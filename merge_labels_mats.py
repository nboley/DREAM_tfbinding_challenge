import os, sys
import gzip

from itertools import izip
from collections import defaultdict

import numpy as np

from grit.lib.multiprocessing_utils import run_in_parallel

regions="/mnt/data/TF_binding/DREAM_challenge/annotations/train_regions.bed.gz"
labels_dir = "/mnt/data/TF_binding/DREAM_challenge/chipseq/peaks/label_mats/"

def build_labels_tsv(sample, factors_and_fnames):
    regions_fp = gzip.open(regions)

    print "processing", sample
    factors_and_fnames.sort()
    all_data = []
    for factor, fname in factors_and_fnames:
        all_data.append(np.load(fname))
        #print sample, factor, all_data[-1].shape
    all_data = np.vstack(all_data).T
    # B, U, A
    with open("%s.labels.tsv" % sample, "w") as ofp:
        # write the header
        ofp.write("\t".join(
            ["chr", "start", "stop"] + [x[0] for x in factors_and_fnames])
                  + "\n")
        for i, (region_line, row) in enumerate(izip(regions_fp, all_data)):
            if i%1000000 == 0:
                print i, all_data.shape[0],
                print region_line.strip() + "\t" + "\t".join("UBA"[x] for x in row)
            ofp.write(region_line.strip() + "\t" + "\t".join("UBA"[x] for x in row) + "\n")
    return 

def main():    
    sample_grpd_labels = defaultdict(list)
    for fname in os.listdir(labels_dir):
        if not fname.endswith(".npy"): continue
        factor, sample_type, _ = fname.split(".")
        sample_grpd_labels[sample_type].append(
            (factor, os.path.join(labels_dir, fname)))
    
    args = list(sample_grpd_labels.iteritems())
    run_in_parallel(16, build_labels_tsv, args)
    return
    
if __name__ == '__main__':
    main()
