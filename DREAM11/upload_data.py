import os

import synapseclient
from synapseclient import Project, Folder, File, Link

syn = synapseclient.Synapse()
syn.login('nboley', 'Sp33d427')

chipseq_bams_id = "syn6153155"
chipseq_conservative_peaks_id = "syn6153159"
chipseq_relaxed_peaks_id = "syn6153158"
chipseq_region_labels_id = "syn6153156"
chipseq_region_scores_id = "syn6153160"

expression_data_id = "syn6153148"

def upload_dnase_idr_peaks():
    dnase_idr_peaks_id = "syn6153193"
    idr_peaks_path = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/idr"
    parent = syn.get(dnase_idr_peaks_id)
    for fname in os.listdir(idr_peaks_path):
        data = File(os.path.join(idr_peaks_path, fname), parent=parent)
        data.annotations['sample_type_id'] = fname.split(".")[0]
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_naive_overlap_peaks():
    dnase_naive_overlap_peaks_id = "syn6153192"
    naive_overlap_peaks_path = "/mnt/data/TF_binding/DREAM_challenge/DNASE/peaks/naive_overlap"
    parent = syn.get(dnase_naive_overlap_peaks_id)
    for fname in os.listdir(naive_overlap_peaks_path):
        data = File(os.path.join(naive_overlap_peaks_path, fname),parent=parent)
        data.annotations['sample_type_id'] = fname.split(".")[0]
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_fold_coverage_wiggles():
    dnase_wigs_id = "syn6153151"
    fold_coverage_wiggles_path = "/mnt/data/TF_binding/DREAM_challenge/DNASE/fold_coverage_wiggles/"
    parent = syn.get(dnase_wigs_id)
    for fname in os.listdir(fold_coverage_wiggles_path):
        data = File(os.path.join(fold_coverage_wiggles_path, fname),parent=parent)
        data.annotations['sample_type_id'] = fname.split(".")[0]
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_bams():
    dnase_bams_id = "syn6153150"
    dnase_bams_path = "/mnt/data/TF_binding/DREAM_challenge/DNASE/bams/"
    parent = syn.get(dnase_bams_id)
    for fname in os.listdir(dnase_bams_path):
        data = File(os.path.join(dnase_bams_path, fname),parent=parent)
        data.annotations['sample_type_id'] = fname.split(".")[0]
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_data():
    upload_dnase_idr_peaks()
    upload_dnase_naive_overlap_peaks()
    upload_dnase_fold_coverage_wiggles()
    upload_dnase_bams()
    
    pass

def main():
    upload_dnase_data()
    pass

if __name__ == '__main__':
    main()
