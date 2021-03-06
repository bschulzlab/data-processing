import pandas as pd

# These are modifiable parameters for the input
peptide_filename = r""
fdr_filename = r""
output_protein_filename = r""
fdr_threshold = 0.01
###############################


class FDRRow:
    def __init__(self, protein, peptide, charge, fdr):
        self.protein = protein
        self.peptide = peptide
        self.charge = charge
        self.fdr = fdr
        self.fdr_pass = False

    def check_threshold(self, threshold, reverse=None):
        check = self.fdr < threshold
        if not reverse:
            self.fdr_pass = check
        else:
            self.fdr_pass = check == True and reverse >= threshold
        return self.fdr_pass

# Loading files into pandas dataframe
fdr_df = pd.read_csv(fdr_filename, sep="\t")

# Load peptide file with values from protein name, peptide sequence and precursor charge as index
peptide_df = pd.read_csv(peptide_filename, sep="\t", index_col=[0, 1, 3])

sample_column_order = peptide_df.columns[2:]

# Assuming that sample column order is consistent
fdr_columns_number = len(fdr_df.columns)

# A temporary dictionary that contain the sum of proteins in each sample
temp_protein_dict = {}

# Iterate through fdr files getting only non-decoy value
for ind, row in fdr_df.iterrows():
    if not row["Decoy"]:
        # If row is not decoy, iterate through each column from column 7 (first sample column)
        for i in range(7, fdr_columns_number, 1):
            current_row = FDRRow(row["Protein"], row["Peptide"], row["Precursor Charge"], row[fdr_df.columns[i]])
            if current_row.check_threshold(fdr_threshold):
                # i-5 is the column number of the current sample
                column_name = peptide_df.columns[i-5]
                if column_name not in temp_protein_dict:
                    temp_protein_dict[column_name] = {}
                if row["Protein"] not in temp_protein_dict[column_name]:
                    temp_protein_dict[column_name][row["Protein"]] = 0
                # Add value from the protein sample to the value of the value of the same protein in the same sample
                if (current_row.protein, current_row.peptide, current_row.charge) in peptide_df.index:
                        sub_df = peptide_df.loc[(current_row.protein, current_row.peptide, current_row.charge)]
                        if type(sub_df) == pd.Series:
                            temp_protein_dict[column_name][row["Protein"]] += sub_df[column_name]
                        else:
                            for i2, r2 in sub_df.iterrows():
                                temp_protein_dict[column_name][row["Protein"]] += r2[column_name]
                else:
                    print((current_row.protein, current_row.peptide, current_row.charge), "does not exist in peptide file")

# Create dataframe from the sample
proteins = pd.DataFrame(temp_protein_dict)
proteins = proteins[sample_column_order]

# Write out new protein file
proteins.to_csv(output_protein_filename, sep="\t", index_label="Protein")
