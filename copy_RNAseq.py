import os, sys

fnames = """
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E003.H1_Cells.ENCODE.CSHL_Gingeras.ENCSR000COU.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E003.H1_Cells.ENCODE.CSHL_Gingeras.ENCSR000COU.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E017.IMR90_Fetal_Lung_Fibroblasts_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CTQ.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E017.IMR90_Fetal_Lung_Fibroblasts_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CTQ.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E066.Liver.ENCODE.CSHL_Gingeras.ENCSR000AFB.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E066.Liver.ENCODE.CSHL_Gingeras.ENCSR000AFB.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E114.A549_EtOH_pt02pct_Lung_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CON.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E114.A549_EtOH_pt02pct_Lung_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CON.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E116.GM12878_Lymphoblastoid_Cells.ENCODE.CSHL_Gingeras.ENCSR000AED.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E116.GM12878_Lymphoblastoid_Cells.ENCODE.CSHL_Gingeras.ENCSR000AED.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E117.HeLa-S3_Cervical_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPR.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E117.HeLa-S3_Cervical_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPR.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E118.HepG2_Hepatocellular_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPE.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E118.HepG2_Hepatocellular_Carcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPE.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E123.K562_Leukemia_Cells.ENCODE.CSHL_Gingeras.ENCSR000AEM.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E123.K562_Leukemia_Cells.ENCODE.CSHL_Gingeras.ENCSR000AEM.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E303.SK-N-SH_Neuroblastoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CTT.RNA-seq.b3.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E303.SK-N-SH_Neuroblastoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CTT.RNA-seq.b4.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E329.MCF-7_Mammary_Gland_Adenocarcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPT.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E329.MCF-7_Mammary_Gland_Adenocarcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR000CPT.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E347.PC-3_Prostate_Adenocarcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR420NLC.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E347.PC-3_Prostate_Adenocarcinoma_Cell_Line.ENCODE.CSHL_Gingeras.ENCSR420NLC.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E364.iPS_PGP1_Cells.CSHL_Gingeras.ENCSR938LSP.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E364.iPS_PGP1_Cells.CSHL_Gingeras.ENCSR938LSP.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E345.Panc1_Epithelioid_Carcinoma_Cell_Line.ENCODE.HAIB_Myers.ENCSR000BYM.RNA-seq.b2.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E345.Panc1_Epithelioid_Carcinoma_Cell_Line.ENCODE.HAIB_Myers.ENCSR000BYM.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E333.HCT116_Colorectal_Carcinoma_Cell_Line.ENCODE.Caltech_Wold.ENCSR000CWM.RNA-seq.b1.t1/Quant.genes.results
/mnt/lab_data/kundaje/users/dskim89/boostCellTypes/data/rna/quant/rsem/encode/E333.HCT116_Colorectal_Carcinoma_Cell_Line.ENCODE.Caltech_Wold.ENCSR000CWM.RNA-seq.b2.t1/Quant.genes.results
""".split()

cell_name_mapping = {
    'iPS': 'induced_pluripotent_stem_cell',
    'H1_Cells': 'H1-hESC'
}

def main():
    for fname in fnames:
        cell_type = fname.split(".")[1].split("_")[0]
        if cell_type in cell_name_mapping:
            cell_type = cell_name_mapping[cell_type]
        print fname
        if ".b1." in fname or '.b3.' in fname:
            bs_id = 1
        elif ".b2." in fname or '.b4.' in fname:
            bs_id = 2
        else:
            assert False
        ofname = "gene_expression.{}.biorep{}.tsv".format(cell_type, bs_id)
        os.system("cp {} {}".format(fname, ofname))

if __name__ == '__main__':
    main()
