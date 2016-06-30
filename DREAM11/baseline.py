import os, sys
from itertools import izip
import gzip

import numpy as np
from scipy.stats.mstats import mquantiles

import pybedtools

import pysam
import pandas as pd

from pyDNAbinding.binding_model import DNASequence
from pyDNAbinding.DB import load_binding_models_from_db

column_names = [
    "chr", "start", "stop",
    "A549", "H1-hESC", "HeLa-S3", "HepG2", "IMR-90", "K562", "MCF-7"
]

def aggregate_region_scores(
        scores, quantile_probs = [0.99, 0.95, 0.90, 0.75, 0.50]):
    rv = [scores.mean()/len(scores), scores.max()]
    rv.extend(mquantiles(scores, prob=quantile_probs))
    return rv

class TrainingData(object):
    def __len__(self):
        return len(self.data)
    
    @property
    def factors(self):
        return self.data.columns[3:].values

    def iter_seqs(self, fasta_fname):
        genome = pysam.FastaFile(fasta_fname)
        return (
            genome.fetch(contig, start, stop+1).upper()
            for contig, start, stop
            in izip(self.data['chr'],
                    self.data['start']-400,
                    self.data['stop']+400)
        )
            
    def __init__(self, fname, max_n_rows=None):
        self.factor = os.path.basename(fname).split('.')[0]
        self.data = pd.read_table(fname, nrows=max_n_rows, names=column_names)
        self.seqs = None
    
    def score_regions(self, fasta_fname):
        all_agg_scores = []
        model = load_binding_models_from_db(tf_names=[self.factor,])[0]
        for i, seq in enumerate(self.iter_seqs(fasta_fname)):
            if i%1000 == 0: print >> sys.stderr, i, len(self)
            all_agg_scores.append(
                aggregate_region_scores(
                    DNASequence(seq).score_binding_sites(model, 'MAX')
                )
            )
        all_agg_scores = np.array(all_agg_scores)
        return all_agg_scores

def load_DNASE_regions_bed():
    try:
        return pybedtools.BedTool("merged_regions.relaxed.bed")
    # pybedtools returns a value error when it can't open a file. If
    # we can't load the file, create and save it
    except ValueError:
        pass
    beds = []
    dnase_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/naive_overlap/"
    for fname in os.listdir(dnase_path):
        beds.append(pybedtools.BedTool(os.path.join(dnase_path, fname)))
    chipseq_peaks_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/idr/"
    for fname in os.listdir(chipseq_peaks_path):
        beds.append(pybedtools.BedTool(os.path.join(chipseq_peaks_path, fname)))
    print "Starting merge"
    regions = beds[0].cat(*beds[1:])
    regions.saveas("merged_regions.relaxed.bed")
    return regions

def load_dnase_filtered_regions():
    # try to load a cached version of this file, and create it
    # if we can't find one
    cached_fname = os.path.basename(sys.argv[1]) + ".relaxed_dnase_filtered_cache"
    try:
        return pybedtools.BedTool(cached_fname)
    except ValueError:
        pass
    dnase_regions_bed = load_DNASE_regions_bed()
    regions_bed = pybedtools.BedTool(sys.argv[1])
    dnase_filtered_regions_bed = regions_bed.intersect(
        dnase_regions_bed, wa=True)
    dnase_filtered_regions_bed.saveas(cached_fname)
    return dnase_filtered_regions_bed
    
def main():
    # find the regions that overlap DNASE peaks in any cell type
    dnase_filtered_regions_bed = load_dnase_filtered_regions()
    print "Finished loading regions."
    # parse the label data
    print dnase_filtered_regions_bed.fn
    data = TrainingData(dnase_filtered_regions_bed.fn)
    print "Finished loading seqs."
    scores = data.score_regions(sys.argv[2])
    np.save("CTCF.v2.motifscores", scores)
    print scores
    
    # build the score traces

    # fit the gradient boosting

    # score the test set

    # 
    return


if __name__ == '__main__':
    main()
