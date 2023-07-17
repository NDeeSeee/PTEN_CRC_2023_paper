import os
import glob
import pandas as pd
import sys
import math
from scipy.stats import binom


### Input parameters
significance = 0.005
sliding_frame_size = 5

### Input data
mut_frequency = pd.read_csv(sys.argv[1], sep="\t")


### Calculation for mid table
mid_table = []

# you can change max number of mutations
for i in range(2,40):
	num_muts_in_codon = i
	num_codons_with_mut = mut_frequency[mut_frequency["Muts"] == num_muts_in_codon].shape[0]
	num_codons_less_mut = mut_frequency[mut_frequency["Muts"] <= num_muts_in_codon].shape[0]
	sum_codons_less_mut = mut_frequency[mut_frequency["Muts"] <= num_muts_in_codon]["Muts"].sum()
	mid_table.append([num_muts_in_codon, num_codons_with_mut, num_codons_less_mut, sum_codons_less_mut])

mid_table = pd.DataFrame(mid_table, columns=["num_muts_in_codon", "num_codons_with_mut", "num_codons_less_mut", "sum_codons_less_mut"])

### Calculation for right table

right_table = {"P-value for exact match":[], "P-value for that many or more":[]}

for i in range(2, 7):
	right_table["shift_" + str(i)] = []

# define cutoff that closest to significance

for i, row in mid_table.iterrows():
	p_val_exact_match = binom.pmf(row["num_muts_in_codon"] - 1, row["sum_codons_less_mut"] - 1, 1/row["num_codons_less_mut"])
	right_table["P-value for exact match"].append(p_val_exact_match)
	p_val_cummulative = binom.cdf(row["sum_codons_less_mut"] - row["num_muts_in_codon"], row["sum_codons_less_mut"] - 1, 1.0 - 1.0/row["num_codons_less_mut"])
	right_table["P-value for that many or more"].append(p_val_cummulative)
	
cutoff = None
for i, pval in enumerate(right_table["P-value for that many or more"]):
	if pval < significance:
		cutoff = i + 2
		break

cutoff = mid_table[mid_table["num_muts_in_codon"] == cutoff].to_dict(orient="list")
cutoff["num_muts_in_codon"][0] -= 1

for i in range(2, 40):
	for j in range(2, 6):
		p_val_cummulative = binom.cdf(cutoff["sum_codons_less_mut"][0] - i, cutoff["sum_codons_less_mut"][0] - 1, 1.0 - float(j)/cutoff["num_codons_less_mut"][0])
		right_table["shift_" + str(j)].append(p_val_cummulative)
	j += 1
	p_val_cummulative = binom.cdf(cutoff["sum_codons_less_mut"][0] - i, cutoff["sum_codons_less_mut"][0], 1.0 - float(j)/cutoff["num_codons_less_mut"][0])
	right_table["shift_" + str(j)].append(p_val_cummulative)
	
right_table = pd.DataFrame(right_table)

lowest_values = {}

for j in range(2,7):
	for i, pval in enumerate(right_table["shift_" + str(j)].to_list()):
		if pval < significance:
			lowest_values[j] = i + 2
			break

lows = [""] + [cutoff["num_muts_in_codon"][0] + 1] + list(lowest_values.values())
lows = [str(elem) for elem in lows]

right_table = right_table.append({col:lows[i] for i, col in enumerate(right_table.columns)}, ignore_index=True)

### Left table
mutations = mut_frequency["Muts"].to_list()
sum_frames = []
count_frames = []

for i in range(0, len(mutations)):
	slice = [val for val in mutations[i:i+sliding_frame_size] if val <= cutoff["num_muts_in_codon"][0]]
	sum_frames.append(sum(slice))
	count_frames.append(len(slice))

mut_frequency["Sum of mut"] = sum_frames
mut_frequency["Count of mut"] = count_frames


vals = []

for i in range(len(sum_frames)):
	try:
		value = right_table.loc[sum_frames[i] - 2, "shift_" + str(count_frames[i])]
		value = -math.log(value, 10)
	except:
		value = right_table.loc[sum_frames[i] - 2, "P-value for that many or more"]
	vals.append(value)


hits = []

for val in vals:
	if val == "None":
		hits.append("")
	elif val >= -math.log(0.005, 10):
		hits.append("Hit")
	else:
		hits.append("")

mut_frequency["Hit or miss"] = hits
mut_frequency["='-lg(p-value)'"] = vals


out = pd.concat([mut_frequency, pd.DataFrame([""], columns=[" "]) , mid_table, pd.DataFrame([""], columns=["  "]), right_table], axis=1)
out = out.fillna("")
out.to_csv(sys.argv[2], sep="\t", index=False)

# mut_frequency.to_csv(sys.argv[2], sep="\t", index=False)

