import pandas as pd

#function to read TCGA data and create mutation matrix. returns the mutation matrix
def get_mutation_matrix_from_maf(maf_file_name):
  maf = pd.read_csv(maf_file_name, sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
  mut = pd.crosstab(maf.Tumor_Sample_Barcode, maf.Hugo_Symbol).clip(upper=1)
  return mut


