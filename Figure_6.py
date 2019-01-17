import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from Bielefeld_dataset_basic_python_dicts import *

def openReadableCSV(filename):#muss für jede Iteration neu kreiert werden
    readableCSV=csv.reader(open(str(filename),newline=''),delimiter='\t')
    return readableCSV
def openWriteableCSV(filename):
    writeableCSV=csv.writer(open(str(filename),'w+',newline=''),delimiter='\t')
    return writeableCSV
    
if not os.path.exists(os.path.join(".", "Figures")):
    os.mkdir(os.path.join(".", "Figures"))

geneDict = get_gene_BS_sum_dict()
#open file and assign gene name in a set and chromosome of not before saved genes in a list
BSFile = openReadableCSV(os.path.join(".","source_files", "all_data.txt"))
usedGeneNames = set()
chromosomeList = []
for row in BSFile:
    if BSFile.line_num != 1:
        if extractEPDGeneName(row[0]) not in usedGeneNames:
            usedGeneNames.add(extractEPDGeneName(row[0]))
            chromosomeList.append(row[3])
# create bar chart for fig 6
fig, ax = plt.subplots()
ax.set(xlabel = "Chromosome", ylabel = "Frequency of Genes")
# set xticks
xticksList = []
for i in range(1, 23):
    xticksList.append("chr"+str(i))
for i in ["chrX", "chrY"]:
    xticksList.append(i)
plt.xticks(np.arange(24), xticksList, rotation = 60, size = "small")
frequencies = []
exp_valFrequencies = []
for i in xticksList:
    frequencies.append(chromosomeList.count(i))
    exp_valFrequencies.append(0)
for gene in geneDict:
    if geneDict[gene]["geneInfo"]["exp_val"] == "Yes":
        exp_valFrequencies[xticksList.index(geneDict[gene]["BS_data"][0]["Chromosome"])] += 1
plt.bar(range(0, 24, 1), frequencies, label = "Amount of NF-$\kappa$B target genes per chromosome")
plt.bar(range(0, 24, 1), exp_valFrequencies, label = "Verified target genes.")
plt.legend()
plt.tight_layout()
plt.show()
if not os.path.exists(os.path.join(".", "Figures")):
    print("creating Folder for histogram")
    os.mkdir(os.path.join(".", "Figures"))
fig.savefig(os.path.join(".", "Figures","figure5.png"), dpi = 600)
print(len(usedGeneNames))
print(sum(frequencies))

# write the used genes into files
BSFile = openReadableCSV(os.path.join(".", "source_files", "all_data.txt"))
geneNamePerChr = {}
allSet = set()
for chr in xticksList:
    geneNamePerChr[chr] = set()
for i, row in enumerate(BSFile):
    if i == 0:
        continue
    geneNamePerChr[row[3]].add(extractEPDGeneName(row[0]))
    allSet.add(extractEPDGeneName(row[0]))
total = 0
for chr in xticksList:
    actualFile = openWriteableCSV(os.path.join(".", "data", "genes_of_" + chr + ".txt"))
    actualFile.writerow(["Genes of " + chr, "Total Number = " + str(len(geneNamePerChr[chr]))])
    total += len(geneNamePerChr[chr])
    for gene in geneNamePerChr[chr]:
        actualFile.writerow([gene])
print(total)
print(len(allSet))

del xticksList[-2]
print(xticksList)
print("starting analysis")
for gene in geneNamePerChr["chrX"]:
    for chr in xticksList:
        if gene in geneNamePerChr[chr]:
            print(chr)
            print(gene)
