import warnings
import sys
import pandas as pd
from collections import Counter

warnings.filterwarnings('ignore')

def main():
	# import data
	full_table = pd.read_csv(sys.argv[1], sep="\t")
	# replace it
	mut_freq = parse_file_to_dict(sys.argv[2])
	# mut_freq = parse_file_to_dict("mut_in_codons_fm.tsv")
	mtest = pd.read_csv(sys.argv[3], sep="\t")

	# replace codon indicies to mutation number in the random table
	# sort by number of mutations (ascending)
	mtest = match_muts(mtest, mut_freq)

	# add number of mutation for each codon
	full_table = set_mutation_codon(full_table, "Codon_1", mut_freq)
	full_table = set_mutation_codon(full_table, "Codon_2", mut_freq)

	# skip linear hotspots
	skip_codons = parse_hotspots(sys.argv[4])
	# first codon
	full_table["needed"] = full_table["Codon_1"].apply(lambda x: False if x in skip_codons else True)
	# second codon
	full_table["needed1"] = full_table["Codon_2"].apply(lambda x: False if x in skip_codons else True)
	
	# filter by second codon
	full_table = full_table[full_table["needed1"] == True]

	# set mutation=0 for codon 1 if it is in skip_codons
	full_table.loc[full_table["needed"] == False, ["Codon_1__muts"]] = 0

	# filter by coeffs
	full_table = full_table[ (full_table["coef_1"] >= 0.8) & (full_table["coef_2"] >= 0.75) ]

	# grouping table
	x = full_table[["Codon_1", "Codon_2", "needed", "Codon_1__muts", "Codon_2__muts"]]
	by_codon1_sum = x.groupby("Codon_1")["Codon_2__muts"].sum()
	by_codon1_max = x.groupby("Codon_1")["Codon_2__muts"].max()
	by_codon1 = pd.DataFrame({"Codon_1":list(by_codon1_max.index), "sum":list(by_codon1_sum), "max":list(by_codon1_max)})

	# merge x and by_codon1
	by_codon1 = pd.merge(by_codon1, x, on="Codon_1")

	# find max in groups of codon1 -> codon2
	by_codon1["Max_in_group"] = by_codon1[["max", "Codon_1__muts"]].apply(lambda x: x[0] if x[0] >= x[1] else x[1], axis=1)
	by_codon1["Sum_in_group"] = by_codon1[["sum", "Codon_1__muts"]].apply(lambda x: x[0] + x[1], axis=1)

	x = by_codon1[["Codon_1", "needed", "Max_in_group", "Sum_in_group"]]

	# calculate size of groups
	count_codon_1 = dict(Counter(list(x["Codon_1"])))
	x["Count_in_group"] = [count_codon_1[elem] + 1 for elem in list(x["Codon_1"])]
	# remove duplicates
	x = x.drop_duplicates()

	# substruct 1 from needed=false
	x.loc[x["needed"] == False, ["Count_in_group"]] -= 1

	# filter by sum in group and max in group
	x1 = x[(x["Sum_in_group"] > 1) & (x["Sum_in_group"] != x["Max_in_group"])]
	sum_val_list = []
	max_val_list = []
	p_val_list = []
	del x1["needed"]

	# comparison with generated table
	for i, row in x1.iterrows():
		mtn3 = mtest.iloc[:,0:row["Count_in_group"]]
		
		condition_1 = mtn3.max(axis=1)<=row["Max_in_group"]
		condition_2 = mtn3.sum(axis=1)>=row["Sum_in_group"]
		max_val = (condition_1).astype(int).sum()
		sum_val = (condition_1 & condition_2).astype(int).sum()
		sum_val_list.append(sum_val)
		max_val_list.append(max_val)
		p_val_list.append(sum_val/max_val)

	x1["max_test"] = max_val_list
	x1["sum_test"] = sum_val_list
	x1["p_test"] = p_val_list

	final = pd.merge(full_table, x1, on="Codon_1").sort_values(by="p_test")
	# final.to_csv("final.tsv", sep="\t", index=False)
	final.to_csv(sys.argv[5], sep="\t", index=False)
	
def parse_file_to_dict(file):
	with open(file, "r") as f:
		f = f.read().split("\n")[1::]

	return {int(line.split("\t")[0]):int(line.split("\t")[1]) for line in f if "\t" in line}


def match_muts(mtest, mut_freq):

	cols = list(mtest.columns.values)

	for col in cols:
		mtest[col] = mtest[col].apply(lambda x: mut_freq[x] if x in mut_freq else 0)

	total_muts = []

	for index, row in mtest.iterrows():
		num_zeros = 0
		for elem in row:
			if elem == 0:
				num_zeros += 1
		total_muts.append(len(cols) - num_zeros)

	mtest["total_muts"] = total_muts
	mtest = mtest.sort_values(by="total_muts")
	del mtest["total_muts"]
	
	return mtest

def set_mutation_codon(df, codon_name, mut_freq):
	
	df[codon_name + "__muts"] = df[codon_name].apply(lambda x: mut_freq[x] if x in mut_freq else 0)
	return df

def parse_hotspots(file_name):
	with open(file_name, "r") as f:
		f = f.read().split("\n")

	for line in f:
		if "Linear hotspots" in line:
			hotspots = line.split("\t")[-1].split(" ")

	hotspots = [int(elem) for elem in hotspots]

	return hotspots

main()