import os, sys
from itertools import izip
import gzip
import pysam

import pandas as pd

from pyDNAbinding.binding_model import FixedLengthDNASequences
from pyDNAbinding.DB import load_binding_models_from_db

def optional_gzip_open(fname, mode='r'):
    if fname.split('.')[-1] in ('gz', 'gzip'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)
    
class TrainingData(object):
    @property
    def factors(self):
        return self.data.columns[3:].values

    def init_seqs(self, fasta_fname):
        genome = pysam.FastaFile(fasta_fname)
        seqs_iter = (
            genome.fetch(contig, start, stop+1).upper()
            for contig, start, stop
            in izip(self.data['chr'],
                    self.data['start']-400,
                    self.data['stop']+400) )
        self.seqs = FixedLengthDNASequences(seqs_iter)
        
    def __init__(self, fname, max_n_rows=None):
        self.factor = os.path.basename(fname).split('.')[0]
        self.data = pd.read_table(fname, nrows=max_n_rows)
        self.seqs = None
    
    def score_regions(self):        
        model = load_binding_models_from_db([self.factor,])[0]
        scores = self.seqs.score_binding_sites(model, 'MAX')
        print scores.max(1)
        
def main():
    # parse the label data
    data = TrainingData(sys.argv[1], 10000)
    data.init_seqs(sys.argv[2])
    data.score_regions()
    print data.seqs
    
    # build the score traces

    # fit the gradient boosting

    # score the test set

    # 
    return


if __name__ == '__main__':
    main()
