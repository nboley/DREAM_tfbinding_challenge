import os, sys
import re
from itertools import izip_longest
import random
import math
from bisect import bisect
from collections import namedtuple
import gzip

import numpy as np
from scipy.stats import itemfreq, norm

from sklearn.metrics import roc_auc_score, f1_score, precision_recall_curve, auc

def optional_gzip_open(fname):
    return gzip.open(fname) if fname.endswith(".gz") else open(fname)  

# which measures to use to evaluate teh overall score
MEASURE_NAMES = ['recall_at_10_fdr', 'recall_at_50_fdr', 'auPRC']
ValidationResults = namedtuple('ValidationResults', MEASURE_NAMES)

def recall_at_fdr(y_true, y_score, fdr_cutoff=0.05):
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    fdr = 1- precision
    cutoff_index = next(i for i, x in enumerate(fdr) if x <= fdr_cutoff)
    return recall[cutoff_index]

ClassificationResultData = namedtuple('ClassificationResult', [
    'is_cross_celltype',
    'sample_type', # should be validation or test
    'train_chromosomes',
    'train_samples', 

    'validation_chromosomes',
    'validation_samples', 

    'auROC', 'auPRC', 'F1', 
    'recall_at_25_fdr', 'recall_at_10_fdr', 'recall_at_05_fdr',
    'num_true_positives', 'num_positives',
    'num_true_negatives', 'num_negatives'])

class ClassificationResult(object):
    _fields = ClassificationResultData._fields

    def __iter__(self):
        return iter(getattr(self, field) for field in self._fields)

    def iter_items(self):
        return zip(self._fields, iter(getattr(self, field) for field in self._fields))
    
    def __init__(self, labels, predicted_labels, predicted_prbs,
                 is_cross_celltype=None, sample_type=None,
                 train_chromosomes=None, train_samples=None,
                 validation_chromosomes=None, validation_samples=None):
        # filter out ambiguous labels
        index = labels > -0.5
        predicted_labels = predicted_labels[index]
        predicted_prbs = predicted_prbs[index]
        labels = labels[index]

        self.is_cross_celltype = is_cross_celltype
        self.sample_type = sample_type

        self.train_chromosomes = train_chromosomes
        self.train_samples = train_samples

        self.validation_chromosomes = validation_chromosomes
        self.validation_samples = validation_samples
        
        positives = np.array(labels == 1)
        self.num_true_positives = (predicted_labels[positives] == 1).sum()
        self.num_positives = positives.sum()
        
        negatives = np.array(labels == 0)        
        self.num_true_negatives = (predicted_labels[negatives] == 0).sum()
        self.num_negatives = negatives.sum()

        if positives.sum() + negatives.sum() < len(labels):
            raise ValueError, "All labels must be either 0 or +1"
        
        try: self.auROC = roc_auc_score(positives, predicted_prbs)
        except ValueError: self.auROC = float('NaN')
        precision, recall, _ = precision_recall_curve(positives, predicted_prbs)
        prc = np.array([recall,precision])
        self.auPRC = auc(recall, precision)
        self.F1 = f1_score(positives, predicted_labels)
        self.recall_at_50_fdr = recall_at_fdr(labels, predicted_prbs, fdr_cutoff=0.50)
        self.recall_at_25_fdr = recall_at_fdr(labels, predicted_prbs, fdr_cutoff=0.25)
        self.recall_at_10_fdr = recall_at_fdr(labels, predicted_prbs, fdr_cutoff=0.10)
        self.recall_at_05_fdr = recall_at_fdr(labels, predicted_prbs, fdr_cutoff=0.05)

        return

    @property
    def positive_accuracy(self):
        return float(self.num_true_positives)/(1e-6 + self.num_positives)

    @property
    def negative_accuracy(self):
        return float(self.num_true_negatives)/(1e-6 + self.num_negatives)

    @property
    def balanced_accuracy(self):
        return (self.positive_accuracy + self.negative_accuracy)/2    

    def iter_numerical_results(self):
        for key, val in self.iter_items():
            try: _ = float(val) 
            except TypeError: continue
            yield key, val
        return

    def __str__(self):
        rv = []
        if self.train_samples is not None:
            rv.append("Train Samples: %s\n" % self.train_samples)
        if self.train_chromosomes is not None:
            rv.append("Train Chromosomes: %s\n" % self.train_chromosomes)
        if self.validation_samples is not None:
            rv.append("Validation Samples: %s\n" % self.validation_samples)
        if self.validation_chromosomes is not None:
            rv.append("Validation Chromosomes: %s\n" % self.validation_chromosomes)
        rv.append("Bal Acc: %.3f" % self.balanced_accuracy )
        rv.append("auROC: %.3f" % self.auROC)
        rv.append("auPRC: %.3f" % self.auPRC)
        rv.append("F1: %.3f" % self.F1)
        rv.append("Re@0.50 FDR: %.3f" % self.recall_at_50_fdr)
        rv.append("Re@0.25 FDR: %.3f" % self.recall_at_25_fdr)
        rv.append("Re@0.10 FDR: %.3f" % self.recall_at_10_fdr)
        rv.append("Re@0.05 FDR: %.3f" % self.recall_at_05_fdr)
        rv.append("Positive Accuracy: %.3f (%i/%i)" % (
            self.positive_accuracy, self.num_true_positives,self.num_positives))
        rv.append("Negative Accuracy: %.3f (%i/%i)" % (
            self.negative_accuracy, self.num_true_negatives, self.num_negatives))
        return "\t".join(rv)


def build_sample_test_file(truth_fname, score_column_index, output_fname):
    ofp = open(output_fname, "w")
    with optional_gzip_open(truth_fname) as fp:
        for i, line in enumerate(fp):
            # skip the header
            if i == 0: continue
            if i%1000000 == 0: print "Finished processing line", i
            data = line.split()
            region = data[:3]
            label = data[score_column_index]
            if label == 'U':
                score = random.uniform(0, 0.8)
            elif label == 'A':
                score = random.uniform(0.5, 0.9)
            elif label == 'B':
                score = random.uniform(0.7, 1.0)
            data = region + [score,]
            ofp.write("{}\t{}\t{}\t{}\n".format(*data))
    ofp.close()
    return 

def verify_file_and_build_scores_array(truth_fname, submitted_fname):
    # make sure the region entries are identical and that
    # there is a float for each entry
    truth_fp_iter = iter(optional_gzip_open(truth_fname))
    submitted_fp_iter = iter(optional_gzip_open(submitted_fname))

    # skip the header
    next(truth_fp_iter)
    
    scores = []
    labels = []
    for i, (t_line, s_line) in enumerate(
            izip_longest(truth_fp_iter, submitted_fp_iter)):
        # make sure the files have the same number of lines
        if t_line is None:
            raise ValueError, "The submitted file has more rows than the reference file"
        if s_line is None:
            raise ValueError, "The reference file has more rows than the reference file"

        # parse the truth line
        t_match = re.findall("(\S+\t\d+\t\d+)\t([UAB])\n", t_line)
        assert len(t_match) == 1, "Line %i in the labels file did not match the expected pattern '(\S+\t\d+\t\d+)\t([UAB])\n'" % i

        # parse the submitted file line, raising an error if it doesn't look
        # like expected
        s_match = re.findall("(\S+\t\d+\t\d+)\t(\S+)\n", s_line)
        if len(s_match) != 1:
            raise ValueError("Line %i in submitted file does not conform to the required pattern: '(\S+\t\d+\t\d+)\t(\S+)\n'" % i)
        if t_match[0][0] != s_match[0][0]:
            raise ValueError("Line %i in submitted file does not match line %i in the reference regions file" % (i, i))

        # parse and validate the score
        try:
            score = float(s_match[0][1])
        except ValueError:
            raise ValueError("The score at line %i in the submitted file can not be interpreted as a float" % i)
        scores.append(score)

        # add the label
        label = t_match[0][1]
        if label == 'A':
            labels.append(-1)
        elif label == 'B':
            labels.append(1)
        elif label == 'U':
            labels.append(0)
        else:
            assert False, "Unrecognized label '%s'" % label
    return np.array(labels), np.array(scores)

def iter_block_samples(labels, scores, num_splits=10):
    block_indices = range(0, len(labels), len(labels)/num_splits)
    slices = [slice(*x) for x in zip(block_indices[:-1], block_indices[1:])]
    # leave one slice out
    for slice_index_to_skip in xrange(len(slices)):
        subsampled_scores = np.concatenate([
            scores[x] for i, x in enumerate(slices)
            if i != slice_index_to_skip
        ], axis=0)
        subsampled_labels = np.concatenate([
            labels[x] for i, x in enumerate(slices)
            if i != slice_index_to_skip
        ], axis=0)
        yield subsampled_labels, subsampled_scores
    return

def empirical_p_value(score, other_scores):
    other_scores.sort()
    rank = bisect(other_scores, score)
    return min(0.5, float(rank)/(len(other_scores)+1))

def calc_score_improved_pvalue(labels, scores, previous_scores):
    # bootstrap the validation scores
    results = []
    for i, (ss_labels, ss_scores) in enumerate(
            iter_block_samples(labels, scores)):
        results.append(
            ClassificationResult(ss_labels, ss_scores.round(), ss_scores)
        )

    p_values = []
    for measure_name in MEASURE_NAMES:
        # find the mean and sd for all of the bootstrapped score estimates
        measures = np.array([getattr(x, measure_name) for x in results])
        mean, std = measures.mean(), measures.std()
        p_val = norm.cdf(getattr(previous_scores, measure_name), mean, std)
        # make the bonferroni correction
        p_values.append(p_val)

    # calculate the probability that all of the measures are better than the
    # previous submission 
    p_values = np.array(p_values)
    p_val = (1 - np.product(1-p_values))
    return float(p_val)

def build_combined_score(results, others_results):
    p_values = []
    for measure_name in MEASURE_NAMES:
        p_value = empirical_p_value(
            getattr(results, measure_name),
            [getattr(x, measure_name) for x in others_results]
        )
        p_values.append(p_value)
    return sum(-math.log10(x) for x in p_values)

def verify_and_score_submission(
        ref_fname, submitted_fname, previous_scores, other_participants_results):
    # build the label and score matrices
    labels, scores = verify_file_and_build_scores_array(
        ref_fname, submitted_fname)

    # decide if the score is better than previous scores
    p_val = calc_score_improved_pvalue(labels, scores, previous_scores)
    full_results = ClassificationResult(labels, scores.round(), scores)

    # compute the measures on the full submitted scores
    results = ValidationResults(*[
        getattr(full_results, attr_name) for attr_name in MEASURE_NAMES])

    # compute the combined score
    combined_score = build_combined_score(results, other_participants_results)
    
    return combined_score, p_val, results

def build_test_data():
    recalls_at_10_fdr = [0.654, 0.652, 0.651, 0.694, 0.624]
    recalls_at_50_fdr = [0.662, 0.662, 0.660, 0.672, 0.663]
    auPRCs = [0.679, 0.703, 0.697, 0.701, 0.702]
    #auROCs = [0.978, 0.979, 0.982, 0.978, 0.979]

    other_participants_results = [
        ValidationResults(*x) for x in zip(
            recalls_at_10_fdr, recalls_at_50_fdr, auPRCs)
    ]
    
    previous_scores = ValidationResults(
        random.choice(recalls_at_10_fdr),
        random.choice(recalls_at_50_fdr),
        random.choice(auPRCs)
    )

    return other_participants_results, previous_scores

def main():
    other_participants_results, previous_scores = build_test_data()
    labels_file = sys.argv[1]
    build_sample_test_file(labels_file, 3, "test.scores")
    submission_file = "test.scores"
    combined_score, p_val, results = verify_and_score_submission(
        labels_file, submission_file, previous_scores, other_participants_results)
    print combined_score, p_val, results

if __name__ == '__main__':
    main()
