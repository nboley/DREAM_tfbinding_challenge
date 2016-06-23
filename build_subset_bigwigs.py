import os, sys
base_path = "/mnt/data/TF_binding/DREAM_challenge/chipseq/fold_change_signal/"

def main():
    for fname in os.listdir(base_path):
        ofname = ".".join(fname.split(".")[:3]) + ".train.fc.signal.bw"
        cmd = "bwtool remove mask -inverse train_regions.blacklisted.merged.bed {} {}".format(
            os.path.join(base_path, fname), ofname)
        print cmd

main()
