from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd

# Path to gostats result file.
gostats_result = r""
p_threshold = 0.05

# This will give both bonferroni and fdr_bh corrected p_value in separate columns
methods = ["bonferroni", "fdr_bh"]


results = []
df = pd.read_csv(gostats_result, sep="\t")
for ont, d in df.groupby("Ontology"):
    d = d.sort_values("Pvalue")
    for m in methods:
        r, p_corrected, als, alb = multipletests(d["Pvalue"].values, p_threshold, is_sorted=True, method="bonferroni")
        d["Pvalue_"+m] = p_corrected
    results.append(d)
results = pd.concat(results)
results.to_csv(gostats_result+"_corrected.txt", sep="\t", index=False)
