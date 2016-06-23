import os, sys
from collections import defaultdict

import numpy as np
from scipy.stats import itemfreq

import pybedtools

from lib import build_label_array_from_bed

def _build_regions_labels_from_beds(
        regions_bed, regions_w_flank_bed, feature_bed):
    # we assume that regions_bed and regions_w_flank_bed are sorted
    feature_bed = feature_bed.sort()
    print "Feature bed", feature_bed.fn
    print "Regions bed", regions_bed.fn
    print "Building flank labels ...",
    flank_labels_bed = regions_w_flank_bed.intersect(
            b=feature_bed, c=True, f=1e-12, F=1e-12, e=True)
    flank_labels = build_label_array_from_bed(flank_labels_bed)
    print "Building core labels ...",
    core_labels_bed = regions_bed.intersect(
            b=feature_bed, c=True, f=0.5, F=0.5, e=True)
    core_labels = build_label_array_from_bed(core_labels_bed)
    core_labels_bed.saveas("tmp.test.bed")
    labels = core_labels.copy()
    labels[core_labels != flank_labels] = -1
    print "FINISHED Building labels"
    print itemfreq(labels)
    return labels

def build_regions_peaks_labels(
        regions_bed,
        regions_w_flank_bed,
        optimal_peaks_fname,
        relaxed_peaks_fname):
    optimal_peaks_bed = pybedtools.BedTool(optimal_peaks_fname).sort().merge()
    relaxed_peaks_bed = pybedtools.BedTool(relaxed_peaks_fname).sort().merge()
    optimal_labels = _build_regions_labels_from_beds(
        regions_bed, regions_w_flank_bed, optimal_peaks_bed)
    print optimal_labels

def find_relaxed_peak_bed_tools(tfname):
    path = "/mnt/data/TF_binding/DREAM_challenge/chipseq/peaks/relaxed_peaks/"
    peak_files = defaultdict(list)
    for fname in os.listdir(path):
        file_tfname, sample = fname.split(".")[1:3]
        if file_tfname != tfname: continue
        peak_files[sample].append(os.path.join(path, fname))

    bed_tools = []
    for sample_type, fnames in peak_files.iteritems():
        print sample_type, fnames
    assert False

def main():
    tf_name = sys.argv[1]
    regions_bed = pybedtools.BedTool(sys.argv[2]) #.sort()
    relaxed_peak_bed_tools = find_relaxed_peak_bed_tools(tf_name)
    
    #build_regions_peaks_labels(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
