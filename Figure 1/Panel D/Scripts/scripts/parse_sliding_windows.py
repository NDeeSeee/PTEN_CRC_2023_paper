import sys
import glob

lst = glob.glob(sys.argv[1])

def get_hits(file):
	with open(file, "r") as f:
		f = f.read().split("\n")

	lst = []

	for i, line in enumerate(f):
		if i == 0:
			continue
		if line == "":
			continue
		if line.split("\t")[4] == "Hit":
			lst.append(int(line.split("\t")[0]))

	return lst

dic = {i:0 for i in range(1,190)}

for file in lst:
	freqs = get_hits(file)
	for elem in freqs:
		dic[elem] += 1

for i in range(1,190):
	print(i, dic[i], sep="\t")