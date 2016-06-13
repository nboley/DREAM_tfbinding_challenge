import synapseclient
from synapseclient import Project, Folder, File, Link

syn = synapseclient.Synapse()
syn.login('nboley', 'Sp33d427')

dnase_bams_id = "syn6153150"
dnase_wigs_id = "syn6153151"
dnase_peaks_id = "syn6153153"

chipseq_bams_id = "syn6153155"
chipseq_conservative_peaks_id = "syn6153159"
chipseq_relaxed_peaks_id = "syn6153158"
chipseq_region_labels_id = "syn6153156"
chipseq_region_scores_id = "syn6153160"

expression_data_id = "syn6153148"

def upload_all_data():
    pass

def main():
    pass

if __name__ == '__main__':
    main()
