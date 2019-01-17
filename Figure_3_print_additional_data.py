#get important data!
import csv
import os

with open(os.path.join(".", "source_files", "all_data.txt"), newline = "") as allFile:
	allFile = csv.reader(allFile, delimiter = "\t")
	allFile = list(allFile)
	print("Number of BS: " + str(len(list(allFile))-1))
	header = allFile[0]
	A_range = ["RELA", 9999,0]
	B_range = ["RELB", 9999,0]
	C_range = ["REL", 9999,0]
	for row in allFile[1:]:
		for i in [A_range, B_range, C_range]:
			if row[header.index("CRE")] == i[0]:
				if int(row[header.index("Score")]) > i[2]:
					i[2] = int(row[header.index("Score")])
				if int(row[header.index("Score")]) <= i[1]:
					i[1] = int(row[header.index("Score")])
	for i in [A_range, B_range, C_range]:
		print(i)