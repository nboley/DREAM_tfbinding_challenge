import os, sys
import re
from itertools import izip_longest
import random
import math

import numpy as np
from scipy.stats import itemfreq

from pyDNAbinding.misc import optional_gzip_open
from pyTFbindtools.cross_validation import ClassificationResult

def build_sample_test_file(truth_fname, score_column_index, output_fname):
    ofp = open(output_fname, "w")
    with optional_gzip_open(truth_fname) as fp:
        for i, line in enumerate(fp):
            # skip the header
            if i == 0: continue
            if i%100000 == 0: print i
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
        assert len(t_match) == 1

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

def iter_block_samples(labels, scores, num_splits=20):
    block_indices = range(0, len(labels), len(labels)/20)
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
    rank = max(i for i, ref_score in enumerate(sorted(other_scores))
               if score > ref_score)
    return max(0.5, float(rank)/len(other_scores+1))

def verify_and_score_submission(ref_fname, submitted_fname):
    # find the existing scores
    recalls_at_10_fdr = [0.654, 0.652, 0.651, 0.694, 0.624]
    recalls_at_50_fdr = [0.662, 0.662, 0.660, 0.672, 0.663]
    auPRCs = [0.679, 0.703, 0.697, 0.701, 0.702]
    
    # build the label and score matrices
    labels, scores = verify_file_and_build_scores_array(
        ref_fname, submitted_fname)
    # bootstrap the validation scores
    results = []
    for i, (ss_labels, ss_scores) in enumerate(
            iter_block_samples(labels, scores)):
        results.append(
            ClassificationResult(ss_labels, ss_scores.round(), ss_scores)
        )
        
    # calculate the overall score
    re_at_fdr_10p_epvs = [
        empirical_p_value(x.recall_at_10_fdr, recalls_at_10_fdr)
        for x in results]
    re_at_fdr_50p_epvs = [
        empirical_p_value(x.recall_at_50_fdr, recalls_at_50_fdr)
        for x in results]
    auPRC_epvs = [empirical_p_value(x.auPRC, auPRCs) for x in results]
    scores = []
    for data in zip(auPRC_epvs, re_at_fdr_50p_epvs, re_at_fdr_10p_epvs):
        scores.append(sum(-math.log(x) for x in data))
    return scores

def main():
    #build_sample_test_file(sys.argv[1], 3, "test.scores")
    res = verify_and_score_submission(sys.argv[1], sys.argv[2])
    print res

main()
