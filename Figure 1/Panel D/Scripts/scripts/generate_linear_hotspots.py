import sys
import collections
from scipy.stats import binom

# Usage:
# python binom_calc.py <mutations frequency by codons> <cutoff (default = 0.005)>

def main():

	# Define constants
	try:
		CUTOFF = float(sys.argv[2])
	except:
		CUTOFF = 0.005

	MIN_SIGNIFICANT_MUT_NUMBER = None

	# Parse mutation file
	mutations = parse_file(sys.argv[1])
	mut_frequency_dist = collections.Counter(mutations.values())

	# List for output
	out = ["Number of mutations in codon	Count of codons with that many mut	Count of codons with less than that many mut	Number of mut in those codons combined	probability for exact match	probability for that many or more"]
		
	# Find first p-value < CUTOFF value
	# i is number of mutations
	# threshold 20 is custom; you can change it
	for i in range(2, 1000):
		# number of successes (i.e. mutation_frequency)
		number_s = i
		# sum of all mutation less or esqual i
		trials = calc_sum(mutations, i)
		# number of codons that have less or esqual mutations to i
		probability_s = calc_count(mutations, i)
		# probability mass function
		pmf = binom.pmf(number_s-1,trials-1,1/probability_s, loc=0)
		# cumulative distribution function
		cdf = binom.cdf(trials-number_s,trials-1,1-1/probability_s, loc=0)
		# prepare line for output
		iter_line = [number_s, trials, probability_s, pmf, cdf]
		# convert to str
		iter_line = [str(elem) for elem in iter_line]
		# append line to final list
		out.append("\t".join(iter_line))
		# Define minimal significant number of mutations
		if MIN_SIGNIFICANT_MUT_NUMBER is None and cdf < CUTOFF:
			MIN_SIGNIFICANT_MUT_NUMBER = number_s
			break

	# Add minimal mut number and linear hotspots
	if MIN_SIGNIFICANT_MUT_NUMBER is None:
		out = ["Minimal mut number\t#NA", "Linear hotspots\t#NA"] + out
	else:
		linear_hotspots = calc_linear_hotspot(mutations, MIN_SIGNIFICANT_MUT_NUMBER)
		linear_hotspots = " ".join([str(elem) for elem in linear_hotspots])
		out = ["Minimal mut number\t" + str(MIN_SIGNIFICANT_MUT_NUMBER), "Linear hotspots\t" + linear_hotspots] + out

	out += ["Threshold\t" + str(CUTOFF)]

	print("\n".join(out))

	
def calc_linear_hotspot(muts, cutoff):
	linear_hotspots = []

	for codon, mut_freq in muts.items():
		if mut_freq >= cutoff:
			linear_hotspots.append(codon)

	return linear_hotspots


def calc_count(muts, cutoff):
	count = 0
	for codon, mut_freq in muts.items():
		if mut_freq <= cutoff:
			count += 1

	return count

def calc_sum(muts, cutoff):
	sum_muts = 0

	for codon, mut_freq in muts.items():
		if mut_freq <= cutoff:
			sum_muts += mut_freq

	return sum_muts


def parse_file(filename):

	# open file
	with open(filename, "r") as f:
		f = f.read().split("\n")[1::]

	muts = {}

	for line in f:
		if line == "":
			continue
		line = line.split("\t")
		muts[int(line[0])] = int(line[-1])

	return muts

main()