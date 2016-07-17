import os, sys
from itertools import izip
import gzip
import tempfile
import hashlib
from subprocess import Popen, PIPE
from collections import namedtuple

import numpy as np
from scipy.stats import itemfreq
from scipy.stats.mstats import mquantiles

import h5py

import pybedtools

import pysam
import pandas as pd

from bw import BigWig

from pyTFbindtools.cross_validation import ClassificationResult

from pyDNAbinding.binding_model import DNASequence, PWMBindingModel
from pyDNAbinding.DB import (
    load_binding_models_from_db, NoBindingModelsFoundError, load_all_pwms_from_db)

GenomicRegion = namedtuple('GenomicRegion', ['contig', 'start', 'stop'])

def load_TAF1_binding_model():
    with open("TAF1_motif.txt") as fp:
        data = []
        for line in fp:
            if line.startswith(">"): continue
            data.append(line.split()[1:])
    data = np.array(data, dtype='float')
    return PWMBindingModel(data)

aggregate_region_scores_labels = [
    "mean", "max", "q99", "q95", "q90", "q75", "q50"]
def aggregate_region_scores(
        scores, quantile_probs = [0.99, 0.95, 0.90, 0.75, 0.50]):
    rv = [scores.mean()/len(scores), scores.max()]
    rv.extend(mquantiles(scores, prob=quantile_probs))
    return rv

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

class LabelData(object):
    def __len__(self):
        return len(self.data)
    
    @property
    def samples(self):
        return self.data.columns[3:].values

    def build_integer_labels(self, label_index):
        labels = np.zeros(self.data.ix[:,3+label_index].shape, dtype=int)
        labels[np.array(self.data.ix[:,3+label_index] == 'A', dtype=bool)] = -1
        labels[np.array(self.data.ix[:,3+label_index] == 'B', dtype=bool)] = 1
        return labels

    def build_train_array(self, normalize=True):
        """Unstack the data frame accessibility score and labels.

        """
        accessibility_score_start_index = 3
        predictors = []
        for i in xrange(accessibility_score_start_index,
                        accessibility_score_start_index+len(self.samples)):
            predictors.append(self.dataframe.ix[:,(0,1,2,i)])
            
        pass
    
    def iter_regions(self, flank_size=0):
        for contig, start, stop in izip(
                self.data['chr'],
                self.data['start']-flank_size,
                self.data['stop']+flank_size):
            yield GenomicRegion(contig, start, stop)
        return
    
    def iter_seqs(self, fasta_fname):
        genome = pysam.FastaFile(fasta_fname)
        return (
            genome.fetch(contig, start, stop+1).upper()
            for contig, start, stop
            in self.iter_regions(flank_size=400)
        )

    def _init_header_data(self, labels_fname):
        with gzip.open(labels_fname) as fp:
            header_line = next(iter(fp))
        header_data = header_line.strip().split("\t")
        if header_data[:3] != ['chr', 'start', 'stop']:
            raise ValueError(
                "Unrecognized header line: '%s'" % header_line.strip())
        self.header = header_data
        return

    def __hash__(self):
        if self._hash is None:
            self._hash = abs(hash(hashlib.md5(str((
                md5(self.labels_fname),
                None if self.regions_fname is None else md5(self.regions_fname),
                self.max_n_rows
            ))).hexdigest()))
        return self._hash

    @property
    def cached_fname(self):
        return "labeldata.%s.%s.obj" % (self.factor, hash(self))

    def _build_dataframe(self):
        # if filter_regions is specified, then restrict the labels to
        # regions that overlap these
        if self.regions_fname is not None:
            # load the data into a bed file
            # zcat {fname} | tail -n +2 | head -n 10000 | \
            #   bedtools intersect -wa -a stdin -b {regions_fname} > {output_fname}
            filtered_regions_fp = tempfile.NamedTemporaryFile("w+")
            p1 = Popen(["zcat", self.labels_fname], stdout=PIPE)
            p2 = Popen(["tail", "-n", "+2",], stdout=PIPE, stdin=p1.stdout)
            # check to see if we should limit the numbere of input rows
            p4_input = None
            # if we want to limit the number of rows, then add a call to head
            if self.max_n_rows is not None:
                p3 = Popen(
                    ["head", "-n", str(self.max_n_rows),], stdout=PIPE, stdin=p2.stdout)
                p4_input = p3.stdout
            else:
                p3 = None
                p4_input = p2.stdout
            p4 = Popen(["bedtools", "intersect",
                        "-wa",
                        "-a", "stdin",
                        "-b", self.regions_fname],
                       stdin=p4_input,
                       stdout=filtered_regions_fp
            )
            p4.wait()
            # Allow p* to receive a SIGPIPE if p(*-1) exits.
            p1.terminate()  
            p2.terminate()
            if p3 is not None: p3.terminate()
            # flush the output file cache, and reset the file pointer
            filtered_regions_fp.flush()
            filtered_regions_fp.seek(0)
            return pd.read_table(
                filtered_regions_fp, nrows=self.max_n_rows, names=self.header)            
        else:
            return pd.read_table(
                self.labels_fname, nrows=self.max_n_rows, names=self.header)
    
    def __init__(self,
                 labels_fname,
                 regions_fname=None,
                 max_n_rows=None,
                 load_cached=True):
        self.labels_fname = labels_fname
        self.regions_fname = regions_fname
        self.max_n_rows = max_n_rows
        self._hash = None
        self.load_cached = load_cached
        # extract the sample names from the header
        assert labels_fname.endswith("labels.tsv.gz"), \
            "Unrecognized labels filename '%s'" % regions_fname
        self._init_header_data(labels_fname)
        # extract the factor from the filename
        self.factor = os.path.basename(labels_fname).split('.')[0]

        # if we want to use a cached version...
        if self.load_cached is True:
            try:
                self.h5store = h5py.File(self.cached_fname)
                self.data = pd.read_hdf(self.cached_fname, 'data')
            except KeyError:
                self.data = self._build_dataframe()
                self.data.to_hdf(self.cached_fname, 'data')
                print self.h5store
        else:
            self.data = self._build_dataframe()
        
        return
    
    def build_motif_scores(self, fasta_fname):
        all_agg_scores = np.zeros(
            (len(self), len(aggregate_region_scores_labels)), dtype=float)
        try:
            models = load_binding_models_from_db(tf_names=[self.factor,])
            assert len(models) == 1, "Multiple binding models found for '{}'".format(self.factor)
        except NoBindingModelsFoundError:
            # if we couldnt find a good motif, just find any motif
            # special case TAF1 because it doesnt exist in CISBP
            if self.factor == 'TAF1':
                models = [load_TAF1_binding_model(),]
            else:
                models = load_all_pwms_from_db(tf_names=self.factor)
        model = models[0]
        for i, seq in enumerate(self.iter_seqs(fasta_fname)):
            if i%10000 == 0: print >> sys.stderr, i, len(self)
            all_agg_scores[i,:] = aggregate_region_scores(
                DNASequence(seq).score_binding_sites(model, 'MAX')
            )        
        all_agg_scores = pd.DataFrame(
            all_agg_scores, columns=aggregate_region_scores_labels)
        return all_agg_scores

    def load_or_build_motif_scores(self, fasta_fname):
        try:
            self.motif_scores = pd.read_hdf(self.cached_fname, 'motif_scores')
        except KeyError:
            self.motif_scores = self.build_motif_scores(fasta_fname)
            self.motif_scores.to_hdf(self.cached_fname, 'motif_scores')
        return self.motif_scores

    def build_dnase_fc_scores(self):
        path="/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/fold_coverage_wiggles/"
        scores = np.zeros((len(self), len(self.samples)), dtype=float)
        for sample_i, sample_name in enumerate(self.samples):
            fname = "DNASE.{}.fc.signal.bigwig".format(sample_name)
            b = BigWig(os.path.join(path, fname))
            for region_i, region in enumerate(self.iter_regions()):
                if region_i%1000000 == 0:
                    print "Sample %i/%i, row %i/%i" % (
                        sample_i+1, len(self.samples), region_i, len(self))
                scores[region_i, sample_i] = b.stats(
                    region.contig, region.start, region.stop, 'max')
            b.close()
        return pd.DataFrame(np.nan_to_num(scores), columns=self.samples)
    
    def load_or_build_dnase_fc_scores(self):
        try:
            self.dnase_fc_scores = pd.read_hdf(self.cached_fname, 'dnase_scores')
        except KeyError:
            self.dnase_fc_scores = self.build_dnase_fc_scores()
            self.dnase_fc_scores.to_hdf(self.cached_fname, 'dnase_scores')
        return self.dnase_fc_scores

    def dataframe(self, fasta_fname, normalize=True):
        # load and normalzie the motif and DNASE scores
        motif_scores = self.load_or_build_motif_scores(fasta_fname)
        if normalize:
            motif_scores = (motif_scores - motif_scores.mean())/motif_scores.std()
        dnase_scores = self.load_or_build_dnase_fc_scores()
        if normalize:
            dnase_scores = (dnase_scores - dnase_scores.mean())/dnase_scores.std()

        # concat the DNASE scores into a single column, and repeat the motif
        # scores so that each DNASE entry has a single motif score entry
        motif_scores = pd.concat(
            [motif_scores]*len(dnase_scores.columns), axis=0, ignore_index=True)
        dnase_scores = pd.concat([
            dnase_scores[column] for column in dnase_scores.columns
        ], axis=0, ignore_index=True)
        # concat the dnase scores and motif scores, and then rename the columns
        rv = pd.concat([motif_scores, dnase_scores], axis=1)
        rv.columns = list(motif_scores.columns) + ['dnase_score',]
        return rv
    
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

def load_dnase_filtered_regions(regions_fname):
    assert regions_fname.endswith(".train.labels.tsv.gz"), \
        "Unrecognized labels filename '%s'" % regions_fname
    print regions_fname
    assert False
    # try to load a cached version of this file, and create it
    # if we can't find one
    cached_fname = os.path.basename(regions_fname) + \
                   ".relaxed_dnase_filtered_cache"
    try:
        return pybedtools.BedTool(cached_fname)
    except ValueError:
        pass
    dnase_regions_bed = load_DNASE_regions_bed()
    regions_bed = pybedtools.BedTool(regions_fname)
    dnase_filtered_regions_bed = regions_bed.intersect(
        dnase_regions_bed, wa=True)
    dnase_filtered_regions_bed.saveas(cached_fname)
    return dnase_filtered_regions_bed


def train_model():
    fasta_fname = sys.argv[2]
    train_data = LabelData(
        sys.argv[1], regions_fname=load_DNASE_regions_bed().fn) #, max_n_rows=1000000)

    labels = []
    for sample_index, sample_name in enumerate(train_data.samples):
        labels.append(train_data.build_integer_labels(sample_index))
    labels = np.concatenate(labels, axis=0)
    
    df = train_data.dataframe(fasta_fname)
    
    # filter out ambiguous labels
    df = df.ix[labels > -0.5]
    labels = labels[labels > -0.5]    
    print itemfreq(labels)

    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression, SGDClassifier
    #from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.cross_validation import KFold
    for train_i, test_i in KFold(labels.shape[0], n_folds=5, shuffle=False):
        train_df = np.nan_to_num(df.ix[train_i,:])
        train_labels = labels[train_i]
        test_df = np.nan_to_num(df.ix[test_i,:])
        test_labels = labels[test_i]
        if False:
            ones_labels = np.nonzero(train_labels == 1)[0]
            zeros_labels = np.nonzero(train_labels == 0)[0]
            # subsample the zero labels to be balanced
            balanced_zero_labels = np.random.choice(
                zeros_labels.ravel(), len(ones_labels), replace=False)
            balanced_indices = np.concatenate(
                [balanced_zero_labels, ones_labels], axis=0)
            train_df = train_df[balanced_indices,:]
            train_labels = train_labels[balanced_indices]
            print itemfreq(train_labels)

        
        #mo = RandomForestClassifier(16, n_jobs=16)
        #mo = GradientBoostingClassifier(
        #    max_depth=6, n_estimators=100) #, learning_rate=1.0)
        #mo = LogisticRegression()
        mo = SGDClassifier(
            loss='log', class_weight='balanced', n_jobs=16)
        mo.fit(train_df, train_labels)

        pred = mo.predict(train_df)
        pred_proba = mo.predict_proba(train_df)[:,1]
        print ClassificationResult(train_labels, pred, pred_proba)

        pred = mo.predict(test_df)
        pred_proba = mo.predict_proba(test_df)[:,1]
        print ClassificationResult(labels[test_i], pred, pred_proba)
        print

def generate_label_data():    
    return

def main():
    #generate_label_data()
    train_model()
    
if __name__ == '__main__':
    main()
