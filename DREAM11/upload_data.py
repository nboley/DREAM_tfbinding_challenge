import os

import synapseclient
from synapseclient import Project, Folder, File, Link

syn = synapseclient.Synapse()
syn.login('nboley', 'Sp33d427')


chipseq_region_labels_id = "syn6153156"
chipseq_region_scores_id = "syn6153160"

################################################################################
#
#

def upload_dnase_bams():
    dnase_bams_id = "syn6176232"
    dnase_bams_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/bams/"
    parent = syn.get(dnase_bams_id)
    for fname in os.listdir(dnase_bams_path):
        data = File(os.path.join(dnase_bams_path, fname), parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data.annotations['biological_replicate'] = int(fname.split(".")[2][6:])
        data.annotations['technical_replicate'] = int(fname.split(".")[3][7:])
        data = syn.store(data)

def upload_dnase_fold_coverage_wiggles():
    dnase_wigs_id = "syn6176233"
    fold_coverage_wiggles_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/fold_coverage_wiggles/"
    parent = syn.get(dnase_wigs_id)
    for fname in os.listdir(fold_coverage_wiggles_path):
        data = File(
            os.path.join(fold_coverage_wiggles_path, fname), parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_idr_peaks():
    dnase_idr_peaks_id = "syn6176235"
    idr_peaks_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/idr/"
    parent = syn.get(dnase_idr_peaks_id)
    for fname in os.listdir(idr_peaks_path):
        data = File(os.path.join(idr_peaks_path, fname), parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_naive_overlap_peaks():
    dnase_naive_overlap_peaks_id = "syn6176236"
    naive_overlap_peaks_path = "/mnt/data/TF_binding/DREAM_challenge/public_data/DNASE/peaks/naive_overlap/"
    parent = syn.get(dnase_naive_overlap_peaks_id)
    for fname in os.listdir(naive_overlap_peaks_path):
        data = File(os.path.join(naive_overlap_peaks_path, fname),parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data = syn.store(data)

def upload_dnase_data():
    upload_dnase_idr_peaks()
    upload_dnase_naive_overlap_peaks()
    upload_dnase_fold_coverage_wiggles()
    upload_dnase_bams()
    return

################################################################################
#
#

def upload_chipseq_conservative_peaks():
    folder_id = "syn6181337"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/idr/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        data = File(os.path.join(path, fname),parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data.annotations['factor'] = fname.split(".")[2]
        data = syn.store(data)

def upload_chipseq_relaxed_peaks():
    folder_id = "syn6181338"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/peaks/naive_overlap/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        data = File(os.path.join(path, fname),parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data.annotations['factor'] = fname.split(".")[2]
        data = syn.store(data)

def upload_chipseq_fold_change_wiggles():
    folder_id = "syn6181334"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/fold_change_signal/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        print fname
        if not fname.endswith(".bw"): continue
        data = File(os.path.join(path, fname),parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data.annotations['factor'] = fname.split(".")[2]
        data = syn.store(data)

def upload_chipseq_labels():
    folder_id = "syn6181335"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/chipseq/labels/tsvs/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        data = File(os.path.join(path, fname),parent=parent)
        data.annotations['factor'] = fname.split(".")[0]
        data = syn.store(data)


def upload_chipseq_data():
    upload_chipseq_labels()
    #upload_chipseq_conservative_peaks()
    #upload_chipseq_relaxed_peaks()
    #upload_chipseq_fold_change_wiggles()
    return

################################################################################

def upload_annotations():
    folder_id = "syn6184307"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/annotations/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        data = File(os.path.join(path, fname),parent=parent)
        data = syn.store(data)
    return

def upload_rnaseq_data():
    folder_id = "syn6176231"
    path = "/mnt/data/TF_binding/DREAM_challenge/public_data/RNAseq/"
    parent = syn.get(folder_id)
    for fname in os.listdir(path):
        data = File(os.path.join(path, fname),parent=parent)
        data.annotations['sample_type'] = fname.split(".")[1]
        data.annotations['biological_replicate'] = int(fname.split(".")[2][6:])
        data = syn.store(data)
    return

def upload_essential_data_tar():
    folder_id = "syn6112317"
    fnames = [
        #"/mnt/data/TF_binding/DREAM_challenge/public_data/training_data.annotations.tar",
        #"/mnt/data/TF_binding/DREAM_challenge/public_data/training_data.ChIPseq.tar",
        #"/mnt/data/TF_binding/DREAM_challenge/public_data/training_data.DNASE_bams.tar",
        "/mnt/data/TF_binding/DREAM_challenge/public_data/training_data.DNASE_wo_bams.tar",
        #"/mnt/data/TF_binding/DREAM_challenge/public_data/training_data.RNAseq.tar"
    ]
    parent = syn.get(folder_id)
    for fname in fnames:
        data = File(fname, parent=parent)
        data = syn.store(data)
    return

def main():
    upload_essential_data_tar()
    #upload_rnaseq_data()
    #upload_annotations()
    #upload_dnase_data()
    #upload_chipseq_data()
    pass

if __name__ == '__main__':
    main()
