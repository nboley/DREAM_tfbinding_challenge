import os, sys
import gzip

from itertools import izip
from collections import defaultdict
from subprocess import check_output

import numpy as np

from grit.lib.multiprocessing_utils import run_in_parallel

#from copy_peaks import find_num_peaks

def find_num_peaks(fname):
    return int(check_output("zcat {} | wc -l".format(fname) , shell=True))

def build_labels_tsv(ofname, samples_and_fnames, regions_fname):
    regions_fp = gzip.open(regions_fname)

    print "processing", ofname
    samples_and_fnames.sort()
    all_data = []
    for sample, fname in samples_and_fnames:
        all_data.append(np.load(fname))
    all_data = np.vstack(all_data).T
    # make sure that the label matrix and regions bed file have the same number
    # of peaks
    assert all_data.shape[0] == find_num_peaks(regions_fname), \
        "Mismatched shapes: %s, %s" % (all_data.shape,  find_num_peaks(regions_fname))
    
    # B, U, A
    with gzip.open(ofname + ".gz", "w") as ofp:
        # write the header
        ofp.write("\t".join(
            ["chr", "start", "stop"] + [x[0] for x in samples_and_fnames])
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
        sample, factor, split_name, _, _ = fname.split(".")
        if factor != 'FOXA2' and factor != 'ATF3': continue
        ofname = "%s.%s.labels.tsv" % (factor, split_name)
        sample_grpd_labels[os.path.join(labels_dir, "..", ofname)].append(
            (sample, os.path.join(labels_dir, fname)))

    args = [list(x) + [regions,] for x in sample_grpd_labels.iteritems()]
    #print args
    run_in_parallel(16, build_labels_tsv, args)
    return
    
if __name__ == '__main__':
    main()
