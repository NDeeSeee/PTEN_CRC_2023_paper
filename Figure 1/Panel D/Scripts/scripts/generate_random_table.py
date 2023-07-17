import os
import random

# Set random seed for reproduce reproducibility
random.seed(432451)

# Set output folder
os.system("mkdir -p data/random_tables")

# Generate 1000 tables with random numbers in range [1,190]
for i in range(0,1000):
	# Header
	lst = ["V1	V2	V3	V4	V5	V6	V7	V8	V9	V10	V11	V12	V13	V14"]
	# 10000 random samples
	for j in range(0,10000):
		# Generate one line of 10000
		# Because the largest group have 14 residues
		# We decided to take same value
		line = random.sample(range(1,190), 14)
		# It need for output
		line = [str(elem) for elem in line]
		# Add line to table
		lst.append("\t".join(line))
	# Merge lines
	lst = "\n".join(lst)
	# Write to file
	with open("data/random_tables/" + str(i+1)+".tsv", "w") as f:
		f.write(lst)
		f.close()