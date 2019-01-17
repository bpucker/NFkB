#RGBA colors: RELAred -> cc17177d, RELBgreen -> 129d137c, RELblue -> 1335a870
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

if not os.path.exists(os.path.join(".", "Figures")):
    os.mkdir(os.path.join(".", "Figures"))

with open(os.path.join(".", "source_files", "all_data.txt"), newline = "") as f:
    biDataset = csv.reader(f, delimiter = "\t")
    for row in biDataset:
        firstRow = row
        break
    motifIndex = firstRow.index("CRE")
    chrIndex = firstRow.index("Chromosome")
    chrList = list(range(1, 23))
    chrList.extend(["X", "Y"])
    chrDict = {}
    for i in chrList:
        chrDict["chr" + str(i)] = {"RELA":0, "RELB":0, "REL":0}
    for row in biDataset:
        chrDict[row[chrIndex]][row[motifIndex]] += 1
a_Data = []
b_Data = []
c_Data = []
for chromosome in chrDict:
    a_Data.append(int(chrDict[chromosome]["RELA"]))
    b_Data.append(int(chrDict[chromosome]["RELB"]))
    c_Data.append(int(chrDict[chromosome]["REL"]))
a_Data = np.array(a_Data)
b_Data = np.array(b_Data)
c_Data = np.array(c_Data)
plt.rc("font", size = 9)
fig, ax = plt.subplots()
ax.set(xlabel = "Chromosome", ylabel = "Frequency of NF-$\kappa$B binding sites")
plt.xticks(np.arange(24), chrDict.keys(), rotation = 60, size = 6.5)
plt.bar(range(24), c_Data, label = "Amount of c-Rel binding sites.", color = "#2035a880", edgecolor = "#1335a870", bottom = 0)
plt.bar(range(24), b_Data, label = "Amount of RelB binding sites.", color = "#129d138c", edgecolor = "#129d138c", bottom = c_Data)
plt.bar(range(24), a_Data, label = "Amount of RelA binding sites.", color = "#cc17177d", edgecolor = "#cc17176d", bottom = b_Data + c_Data)
plt.legend()
fig.savefig(os.path.join(".", "Figures", "Figure_3a.png"), dpi = 600)