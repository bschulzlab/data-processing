import pandas as pd
import numpy as np

# Replacing peptide intensity values of peptide ions above the fdr threshold with a blank cell
# Input file path for peptide and fdr tabulated text files. For the script to run correctly, sample order has to be the same in both files.
# fdr_threshold variable hold the value in which if fdr is greater, the cell value would be replace with blank
input_peptides = r"C:\Users\localadmin\PycharmProjects\fdr_filter\20190605_Combined_Byonic_Grape_Yeast_ProteinPilot_SixTransitions_Peptides.txt"
input_fdr = r"C:\Users\localadmin\PycharmProjects\fdr_filter\20190605_Combined_Byonic_Grape_Yeast_ProteinPilot_SixTransitions_FDR.txt"
fdr_threshold = 0.01

output_peptides = r"C:\Users\localadmin\PycharmProjects\fdr_filter\20190605_Combined_Byonic_Grape_Yeast_ProteinPilot_SixTransitions_Peptides_output.txt"

# Read both files into separate dataframes and index them by Protein, Peptide and Precursor Charge
df_peptides = pd.read_csv(input_peptides, sep="\t", index_col=[0, 1, 3])
df_fdr = pd.read_csv(input_fdr, sep="\t", index_col=[0, 1, 4])
peptides_columns = len(df_peptides.columns)
fdr_columns = len(df_fdr.columns)

# Iterate through a selection of fdr file where Decoy value is False
for i, r in df_fdr[df_fdr["Decoy"] == False].iterrows():
    # Iterate through columns from where the sample value starts in both dataframes (column 3 for peptide, column 5 for fdr)
    for n in range(4, fdr_columns):
        # If fdr value is greater than threshold, replace with np.nan (blank)
        if r[df_fdr.columns[n]] > 0.01:
            df_peptides.at[i, df_peptides.columns[n - 2]] = np.nan

# Write out to tabulated format.
df_peptides.to_csv(output_peptides, sep="\t")

