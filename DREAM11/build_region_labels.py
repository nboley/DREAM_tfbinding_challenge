import os, sys
from collections import defaultdict
from itertools import izip

import numpy as np
from scipy.stats import itemfreq

import pybedtools

from lib import build_label_array_from_bed

def build_regions_labels_from_beds(
        regions_bed, optimal_peaks_bed, relaxed_peaks_bed):
    # we assume that regions_bed and regions_w_flank_bed are sorted
    #print "Building optimal labels ..."
    optimal_peaks_bed = optimal_peaks_bed.sort()
    optimal_labels_bed = regions_bed.intersect(
            b=optimal_peaks_bed, c=True, f=0.5, F=0.5, e=True)
    optimal_labels = build_label_array_from_bed(optimal_labels_bed)

    #print "Building relaxed labels ..."
    relaxed_peaks_bed = relaxed_peaks_bed.sort()
    relaxed_labels_bed = regions_bed.intersect(
            b=relaxed_peaks_bed, c=True, f=1e-6, F=1e-6, e=True)
    relaxed_labels = build_label_array_from_bed(relaxed_labels_bed)

    labels = optimal_labels.copy()
    labels[optimal_labels != relaxed_labels] = -1
    print "FINISHED Building Merged Labels"
    print itemfreq(labels)
    return labels

def iter_peak_bed_tools(tfname, path):
    peak_files = defaultdict(list)
    for fname in os.listdir(path):
        if not fname.startswith('ChIPseq'): continue
        sample, file_tfname = fname.split(".")[1:3]
        #print file_tfname, tfname, sample
        if file_tfname != tfname: continue
        peak_files[sample].append(os.path.join(path, fname))

    bed_tools = []
    for sample_type, fnames in sorted(peak_files.iteritems()):
        print sample_type, fnames
        if len(fnames) == 1:
            bed_tools.append(
                (sample_type,
                 pybedtools.BedTool(fnames[0]).sort().merge())
            )
        else:
            bed_tools.append(
                (sample_type,
                 pybedtools.BedTool(fnames[0]).cat(*fnames[1:]).sort().merge())
            )
        yield bed_tools[-1]    
    return

def old_main():
    tf_name = sys.argv[1]
    #regions_bed = pybedtools.BedTool(sys.argv[2]).sort()
    #relaxed_peak_bed_tools = iter_peak_bed_tools(
    #    tf_name, "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/naive_overlap/")
    #optimal_peak_bed_tools = iter_peak_bed_tools(
    #    tf_name, "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/idr/")

    relaxed_peak_bed_tools = iter_peak_bed_tools(
        tf_name, "/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/naive_overlap/")
    optimal_peak_bed_tools = iter_peak_bed_tools(
        tf_name, "/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/idr/")
    
    for (sample_type_1, relaxed_bed_tool), (sample_type_2, optimal_bed_tool) in izip(
            relaxed_peak_bed_tools, optimal_peak_bed_tools):
        print sample_type_1, sample_type_2
        continue
        print "Processing ", sample_type_1, sample_type_2
        assert sample_type_1 == sample_type_2
        labels = build_regions_labels_from_beds(
            regions_bed, relaxed_bed_tool, optimal_bed_tool)
        ofname = "%s.%s" % (tf_name, sample_type_1)
        print "Saving to %s" % ofname
        np.save(ofname, labels)

def main():
    tf_name = sys.argv[1]
    sample = sys.argv[2]
    regions_bed = pybedtools.BedTool(sys.argv[3])#.sort()

    relaxed_pk_fname = "/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/naive_overlap/ChIPseq.{sample}.{tf_name}.relaxed.narrowPeak.gz".format(tf_name=tf_name, sample=sample)
    relaxed_peak_bed_tool = pybedtools.BedTool(relaxed_pk_fname)
    print relaxed_peak_bed_tool

    optimal_pk_fname = "/mnt/data/TF_binding/DREAM_challenge/private_data/chipseq/peaks/idr/ChIPseq.{sample}.{tf_name}.conservative.narrowPeak.gz".format(tf_name=tf_name, sample=sample)
    optimal_peak_bed_tool = pybedtools.BedTool(relaxed_pk_fname)
    print optimal_peak_bed_tool
    return
    print "Processing ", tf_name, sample
    labels = build_regions_labels_from_beds(
        regions_bed, relaxed_bed_tool, optimal_bed_tool)
    ofname = "%s.%s" % (tf_name, sample_type_1)
    print "Saving to %s" % ofname
    np.save(ofname, labels)

    pass

if __name__ == '__main__':
    main()
