import sys
import glob

def parse_file(filename):
	dict_val = {}

	with open(filename, "r") as f:
		f = f.read().split("\n")

	for i in range(1, len(f)):
		if f[i] == "":
			continue
		line = f[i].split("\t")
		codon_1 = line[0]
		codon_2 = line[2]
		p_val = line[-1]
		dict_val[codon_1+"_"+codon_2] = p_val

	return dict_val


results = {}
lst_out = []

i = 0

for file in glob.glob(sys.argv[1]):
	sample = parse_file(file)
	results[int(file.split("/")[-1].split(".")[0])] = sample
	lst_out.append(list(sample.keys()))
	if int(file.split("/")[-1].split(".")[0]) > i:
		i = int(file.split("/")[-1].split(".")[0])

keys = set([key for lst_in in lst_out for key in lst_in])


header = "codon1\tcodon2"

for i in range(0, i):
	header += "\t" + str(i+1)

print(header)

for key in keys:
	line = key.replace("_","\t")
	for j in range(0, i+1):
		try:
			line += "\t" + str(results[j+1][key])
		except:
			line += "\tNaN"
	print(line)

