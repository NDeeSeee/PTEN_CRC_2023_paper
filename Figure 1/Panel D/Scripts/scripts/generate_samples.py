import os
import sys
import random

def parse_file(filename):
	with open(filename, "r") as f:
		f = f.read().split("\n")

	dic = {}

	for i, line in enumerate(f):
		if i == 0:
			continue
		line = line.split("\t")
		if line[0] in dic:
			dic[line[0]].append(line[17])
		else:
			dic[line[0]] = [line[17]]

	lst = [value for key, value in dic.items()]

	return lst


def main():

	random.seed(65748)

	# open fmi mutations file
	# mutations =  read_tsv("mutations.tsv")
	# mutations =  parse_file()
	# create sample of each mutations
	samples = parse_file(sys.argv[1])

	# print(samples)

	# 1/3 of all mutations
	sample_size = len(samples) // 3
	# number of repeats
	sample_time = 1000
	# get small set
	resamples = resampling(samples, sample_size, sample_time)

	os.system("rm -rf data/resampling/")
	os.system("mkdir -p data/resampling/")
	for i, sample in enumerate(resamples):
		sample = calculate_number(sample)
		with open("data/resampling/" + str(i+1) + ".tsv", "w") as f:
			string = "Codon	Muts\n"
			for res_id, muts in sample.items():
				string += res_id + "\t" + str(muts) + "\n"
			f.write(string)
			f.close()


def read_tsv(filename):
	mut_dict = {}
	with open(filename, "r") as f:
		f = f.read().split("\n")

		for i in range(1,len(f)):
			line = f[i]
			mut, count = line.split("\t")
			mut_dict[mut] = int(count)

		return mut_dict

def create_samples(mutations):
	mut_list = []

	for mut_name, mut_number in mutations.items():
		for i in range(mut_number):
			mut_list.append(mut_name)

	return mut_list

def resampling(samples, size, time):

	random.shuffle(samples)
	resamples = []

	for i in range(time):
		random_samples = random.sample(samples, size)
		resamples.append(random_samples)

	return resamples

def calculate_number(sample):
	counter = {}

	for i in range(1,190):
		counter[str(i)] = 0

	for lst in sample:
		for elem in lst:
			counter[elem] += 1

	return dict(counter)

main()