from Bielefeld_dataset_basic_python_dicts import *
import matplotlib.pyplot as plt
import numpy as np
import os
import csv


if not os.path.exists(os.path.join(".", "Figures")):
    os.mkdir(os.path.join(".", "Figures"))

# collect data
promDict = get_prom_BS_sum_dict()
# colorDict contains tuples with (color, edgecolor)
colorDict = {"RELA":("#cc17177d", "#cc17176d"), 
            "RELB":("#129d138c", "#129d138c"), 
            "c-Rel":("#2035a880", "#1335a870"),
            "all": ("#99999980", "#99999980")}
chrList = list(range(1, 23))
chrList.extend(["X", "Y"])
chrDict = {}
avgScoreDict = {"RELA":[], "RELB":[], "REL":[], "all": []}
for i in chrList:
    chrDict["chr" + str(i)] = {"RELA":[], "RELB":[], "REL":[], "all": []}
for promoter in promDict:
    currentPromoterScores = {"RELA":[], "RELB":[], "REL":[], "all": []}
    for BS in promDict[promoter]["BS_data"]:
        currentPromoterScores["all"].append(BS["Score"])
        avgScoreDict["all"].append(BS["Score"])
        currentPromoterScores[BS["CRE"]].append(BS["Score"])
        avgScoreDict[BS["CRE"]].append(BS["Score"])
    for keys in currentPromoterScores:
        if currentPromoterScores[keys] != []:
            chrDict[BS["Chromosome"]][keys].append(np.mean(currentPromoterScores[keys]))
            chrDict[BS["Chromosome"]][keys].append(np.mean(currentPromoterScores[keys]))
        else:
            chrDict[BS["Chromosome"]][keys].append(0)
del promDict
chrAvgScoreDict = {"RELA":[], "RELB":[], "REL":[], "all": []} # holds chromosomes avg Scores sorted by motif type
for chromosome in chrDict:
    for motifType in chrDict[chromosome]:
        chrDict[chromosome][motifType] = np.mean(chrDict[chromosome][motifType])
        chrAvgScoreDict[motifType].append(chrDict[chromosome][motifType])
for motifType in avgScoreDict:
    avgScoreDict[motifType] = np.mean(avgScoreDict[motifType])
chrAvgScoreDict["c-Rel"] = chrAvgScoreDict.pop("REL")
avgScoreDict["c-Rel"] = avgScoreDict.pop("REL")
avgScoreDict["all"] = avgScoreDict.pop("all")
# create barplots
xticks = ["RELA", "RELB", "c-Rel", "all"]
motifToLegendDict = {"RELA":"RelA", "RELB":"RelB", "c-Rel":"c-Rel", "all":"all"}
plt.figure(figsize = (3, 6))
plt.xlabel("Motif Type", size = 11)
plt.ylabel("Average JASPAR Score of NF-$\kappa$B binding sites", size = 11)
indexes = [1, 1.5, 2, 2.5]
plt.xticks(indexes, [motifToLegendDict[i] for i in xticks], rotation = None, size = 9.5)
plt.yticks(np.arange(0, 401, 50), np.arange(0, 401, 50), size = 9.5)
for index, motifType in enumerate(xticks):
    barHeight = avgScoreDict[motifType]
    plt.bar(indexes[index], barHeight, color = colorDict[motifType][0], 
    edgecolor = colorDict[motifType][1], width = 0.4, label = "Average Score of " + motifToLegendDict[motifType] + " binding sites.")
# plt.legend(loc = "lower center", prop = {"size": 7})
plt.tight_layout()
#plt.show()
plt.savefig(os.path.join(".", "Figures", "BT_figure_3b.png"), dpi = 400)