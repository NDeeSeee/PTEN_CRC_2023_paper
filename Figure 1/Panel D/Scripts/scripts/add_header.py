import glob

def add(file):
	with open(file, "r") as f:
		data = f.read()

	with open(file, "w") as f:
		f.write("Codon	Muts\n" + data)
		f.close()

for file in glob.glob("data/samples/*tsv"):
	add(file)