import glob
import sys

def parse_hotspots(file_name):
	with open(file_name, "r") as f:
		f = f.read().split("\n")

	for line in f:
		if "Linear hotspots" in line:
			hotspots = line.split("\t")[-1].split(" ")

	hotspots = [int(elem) for elem in hotspots]

	return hotspots

dic = {}

for file in glob.glob(sys.argv[1]):
	lins = parse_hotspots(file)
	for elem in lins:
		if elem not in dic:
			dic[elem] = 1
		else:
			dic[elem] += 1


for key, value in dic.items():
	print(key, value, sep="\t")