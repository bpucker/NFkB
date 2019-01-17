import csv
import os
import re
import requests
import sys
from numpy import mean as mean

# define helper functions
def openReadableCSV(filePath):
    readableCSV=csv.reader(open(str(filePath),newline=''),delimiter='\t')
    return readableCSV

def openWriteableCSV(filePath):
    writeableCSV=csv.writer(open(str(filePath),'w+',newline=''),delimiter='\t')
    return writeableCSV
    
def extractEPDGeneName(PromoterID):
    EPDName=re.compile(r'([a-zA-Z0-9@\-.]+)_\d+').findall(PromoterID)
    if EPDName == []:
        print("No EPD Name found in " + PromoterID)
        return "Not found"
    else:
        return EPDName[0]

def getEnsemblFromGeneSymbol(GeneSymbol):
    if isinstance(GeneSymbol, str):
        server = "https://rest.ensembl.org"
        ext = "/xrefs/symbol/homo_sapiens/"+GeneSymbol+"?"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        decoded = r.json()
        if decoded != []:
            EnsemblIDList = []
            for i in range(len(decoded)):
                EnsemblIDList.append(decoded[0]["id"])
        else:
            EnsemblIDList = []
        return EnsemblIDList
    else:
        print("pls dont put in more than one gs at once")
        input()

def enrich_basic_data():
    if not os.path.exists(os.path.join(".", "source_files", "promDataForWang.csv")):
        # create dictionary: {promID : {"Start" : , "Sequence" : }}
        promDataForWang = {}
        promDataFile = openReadableCSV(os.path.join(".", "source_files", "promoters.gff3"))
        for row in promDataFile:
            if promDataFile.line_num < 3:
                continue
            elif promDataFile.line_num == 3:
                startIndex = row.index("Start")
                promIndex = row.index("PID")
            else:
                promDataForWang[row[promIndex]] = {"Start" : row[startIndex]}
        notThereList = []
        with open(os.path.join(".", "source_files", "fasta_sequences_from_-5000_to_+100.fa"), "r") as sequenceFile:
            for line_num, row in enumerate(sequenceFile):
                if row[0] == ">" and line_num != 0:
                    try:
                        promDataForWang[currentPromoter]["Sequence"] = currentSequence
                    except KeyError:
                        notThereList.append(currentPromoter)
                    finally:
                        currentPromoter = re.compile(r'\D\D\d+\s[a-zA-Z0-9@\.-]+_\d+').findall(row)
                        currentPromoter = currentPromoter[0].replace(" ", "_%_")
                        currentSequence = ""
                elif row[0] == ">":
                    currentPromoter = re.compile(r'\D\D\d+\s[a-zA-Z0-9@\.-]+_\d+').findall(row)
                    currentPromoter = currentPromoter[0].replace(" ", "_%_")
                    currentSequence = ""
                elif currentSequence == "":
                    currentSequence = row.replace("\n", "")
                else:
                    currentSequence += row.replace("\n", "")
            try:
                promDataForWang[currentPromoter]["Sequence"] = currentSequence
            except KeyError:
                notThereList.append(currentPromoter)
            sequenceFile.close()
        print(str(len(notThereList)), " promoters where not found.")
        promDataForWangFile = openWriteableCSV(os.path.join(".", "source_files", "promDataForWang.csv"))
        for key in promDataForWang:
            try:
                promDataForWangFile.writerow([key, promDataForWang[key]["Start"], promDataForWang[key]["Sequence"]])
            except Exception as exp:
                print(exp)
                print(key)
                print(promDataForWang[key])
                input()
    wangFile = openReadableCSV(os.path.join(".", "source_files", "promDataForWang.csv"))
    wangDict = {} # structure: {promID : {"Start" : , "Sequence" : }}
    for row in wangFile:
        wangDict[row[0]] = {"Start": int(row[1]), "Sequence": row[2]}
    allDataFile = openReadableCSV(os.path.join(".", "source_files", "all_data.txt"))
    expandedAllDataFile = openWriteableCSV(os.path.join(".", "source_files", "all_data_expanded.txt"))
    iDict = {}
    for row in allDataFile:
        if allDataFile.line_num == 1:
            row.append("BS_Sequence")
            row.append("RELA_central_BS_W")
            row.append("RELA_cW_sur_fit")
            for index, column in enumerate(row):
                iDict[column] = index
        else:
            length = int(row[iDict["End"]]) - int(row[iDict["Start"]])
            readingStartPoint = int(row[iDict["Start"]]) - wangDict[row[iDict["PromoterID"]]]["Start"]
            readingEndPoint = readingStartPoint + length
            BSSequence = wangDict[row[iDict["PromoterID"]]]["Sequence"][readingStartPoint : readingEndPoint]
            row.append(BSSequence)
            if BSSequence[5].upper() in ["A", "T"]:
                row.append(1)
                if BSSequence[3:8].upper() in ["TTTAA", "AAATT", "AATTT", "TTAAA"]:
                    row.append(1)
                else:
                    row.append(0)
            else:
                row.append(0)
                row.append(0)
        expandedAllDataFile.writerow(row)
    
    
def create_foundation_file():
    if (not os.path.exists(os.path.join(".", "source_files","BS_sum_complData.csv"))):
        print("creating foundation file for dictionary creation")
        if( not os.path.exists(os.path.join(".", "source_files","all_data_expanded.txt"))):
            enrich_basic_data()
        # expanded result of the analysis by bpucker as assembled with the scripts on https://github.com/bpucker/NFkB
        BSFile = openReadableCSV(os.path.join(".","source_files", "all_data_expanded.txt"))
        # additional information about the chromosomes provided by bpucker
        promoterFile = openReadableCSV(os.path.join(".","source_files", "promoters.gff3"))
        # mapping files from the EPD FTP server: ftp://ccg.vital-it.ch/epdnew/H_sapiens/006/db/ 
        prom_ens = openReadableCSV(os.path.join(".","source_files", "promoter_ensembl.txt")) 
        cross_ref = openReadableCSV(os.path.join(".","source_files", "cross_references.txt"))
        # target file which serves as the basis of all created python dictionaries
        BS_sumFile = openWriteableCSV(os.path.join("source_files","BS_sum_complData.csv"))
        
        #save data to add in dictionaries
        ensemblDataDict = {} #saves data from cross_ref with ensembl as key and [ens, geneName, refseq, description, exp_val]
        geneNameDataDict = {} #same but with geneName as key

        with open(os.path.join(".", "source_files", "bi_all_overlap_ensembls.csv")) as overlapEnsemblFile:
            overlapEnsembls = set()
            for ensembl in overlapEnsemblFile:
                overlapEnsembls.add(ensembl.replace("\n", ""))
            print("Overlap Ensembls gotten")

        for row in cross_ref:
            ensemblDataDict[row[0]] = row
            geneNameDataDict[row[1]] = row
        promoterDict = {} #gathers all additional (not in bs file) data for every promoter
        for row in prom_ens:
            promoterDict[row[0]] = ensemblDataDict[row[1]]
        del ensemblDataDict
        promEndDict = {}
        #save end of promoter sequence
        for row in promoterFile: 
            if promoterFile.line_num > 3:
                promEndDict[row[7].split("_%_")[1]] = int(row[4]) #end

        dictToWrite={}
        for bs in BSFile: #row in bsfile
            #write header
            if BSFile.line_num == 1:
                col_index={}
                for index,category in enumerate(bs):
                    col_index[category] = index
                col_index["dist_TSS"] = len(col_index.keys())
                BS_sumFile.writerow(list(col_index.keys())[1:])
            #save rows to write in dictToWrite as promoter ID:[bsFileRows, additional Info] 
            else:
                #new entry in dict if promoter ID is not yet in it
                prom_ID = bs[col_index["PromoterID"]].split("_%_")[1]
                #add distance to TSS
                bs.append(str(promEndDict[prom_ID]-int(bs[col_index["End"]])))
                if prom_ID not in dictToWrite:
                    if prom_ID in promoterDict:
                        dictToWrite[prom_ID] = [[bs[1:]], promoterDict[prom_ID]]
                    elif prom_ID[:-2] in geneNameDataDict:
                        dictToWrite[prom_ID] = [[bs[1:]], geneNameDataDict[extractEPDGeneName(prom_ID)]]
                    else:
                        currentEnsemblID = getEnsemblFromGeneSymbol(extractEPDGeneName(prom_ID))
                        if currentEnsemblID != []:
                            dictToWrite[prom_ID] = [[bs[1:]], 
                            ["_%_".join(currentEnsemblID), extractEPDGeneName(prom_ID), "XXX", "XXX"]]
                        elif currentEnsemblID == []:
                            dictToWrite[prom_ID] = [[bs[1:]], ["XXX", extractEPDGeneName(prom_ID), "XXX", "XXX"]]
                        else:
                            print("something is wrong with your code")
                            input()
                else:
                    if len(dictToWrite[prom_ID][0]) > 0:
                        #start heap sort
                        for i, entry in enumerate(dictToWrite[prom_ID][0]):
                            if bs[col_index["Score"]] >= entry[col_index["Score"]-1]:
                                dictToWrite[prom_ID][0].insert(i, bs[1:])
                                break
                            elif i == len(dictToWrite[prom_ID][0])-1:
                                dictToWrite[prom_ID][0].append(bs[1:])
                                break
                    else:
                        print("look at your code")
                        input()
        for promoter in dictToWrite:
            if dictToWrite[promoter][1][0] in overlapEnsembls:
                dictToWrite[promoter][1].append("Yes")
            else:
                dictToWrite[promoter][1].append("No")
            listToWrite = []
            BS_sumFile.writerow([promoter])
            BS_sumFile.writerow(dictToWrite[promoter][1])
            for i in dictToWrite[promoter][0]:
                bsInfostring = "_%_".join(i)
                listToWrite.append(bsInfostring)
            BS_sumFile.writerow(listToWrite)
    else:
        pass
        
        


def get_prom_BS_sum_dict(threshold = 0):
    #returns data as dictionary with promoter ID as key, subdictionary holds the keys geneInfo and BS_data
    #The BS_data key holds a list with one dictionary for every bs accessible over their category
    #The categories are the following: GeneID CRE Chromosome Start End  Score Orientation dist_TSS RELA_central_BS_W ( -> QS_Wang_1)    RELA_cW_sur_fit ( -> QS_Wang_5)
    #The geneInfo key holds a dict with the categories: "ENSEMBL_ID", "GeneName", "RefSeq", "Gene_Description", "exp_val"
    #keys of geneInfo and BS_data give access to string values
    
    #first create foundational file if it doesnt already exitsts in the folder "source_files"
    #Attention! If file is corrupted it needs to be deleted, before repeating this program!
    create_foundation_file()
    
    #create promoter Dictionary based on the foundation file
    BS_sum_file = openReadableCSV(os.path.join(".", "source_files","BS_sum_complData.csv"))
    prom_BS_sum_dict = {}
    for row in BS_sum_file:
        if BS_sum_file.line_num == 1:
            header = row
        elif (BS_sum_file.line_num + 1) % 3 == 0 :
            currentPromID = row[0]
            prom_BS_sum_dict[currentPromID] = {}
        elif (BS_sum_file.line_num % 3) == 0: 
            prom_BS_sum_dict[currentPromID]["geneInfo"] = {}
            for i, category in enumerate(["ENSEMBL_ID", "GeneName", "RefSeq", "Gene_Description", "exp_val"]):
                if category == "ENSEMBL_ID":
                    tmpSet = set()
                    tmpSet.update(row[i].split("_%_"))
                    prom_BS_sum_dict[currentPromID]["geneInfo"][category] = "_%_".join(tmpSet)
                else:
                    try:
                        prom_BS_sum_dict[currentPromID]["geneInfo"][category] = row[i]
                    except:
                        print(i)
                        print(category)
                        print(row)
                        input()
        elif (BS_sum_file.line_num + 2) % 3 == 0 :
            prom_BS_sum_dict[currentPromID]["BS_data"] = []
            for i,content in enumerate(row):
                currentBS = content.split("_%_")
                if int(currentBS[header.index("Score")]) > threshold:
                    prom_BS_sum_dict[currentPromID]["BS_data"].append({})
                    for j in range(len(header)):
                        prom_BS_sum_dict[currentPromID]["BS_data"][i][header[j]] = currentBS[j]
    toConvertToInt = ["Start", "End", "Score", "dist_TSS", "RELA_central_BS_W", "RELA_cW_sur_fit"]
    for promoter in prom_BS_sum_dict:
        for BS in prom_BS_sum_dict[promoter]["BS_data"]:
            for category in toConvertToInt:
                BS[category] = int(BS[category])
    return prom_BS_sum_dict

def get_gene_BS_sum_dict(threshold = 0):
    #creates a dict with genenames as keys, for another dict with the keys "geneInfo" and "BS_data"
    #the only difference to the get_prom_BS_sum_dict function is the added key "num_of_promoters"
    # in the "geneInfo" category
    #The BS_datas get fused without sorting out BSs which occur more than one time in the same gene.
    promDict = get_prom_BS_sum_dict(threshold)
    geneDict = {}
    for promoter in promDict:
        if promDict[promoter]["geneInfo"]["GeneName"] in geneDict:
            #actualize the promoter count
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]]["geneInfo"]["num_of_promoters"] += 1
            #append BSs to "BS_data" dictionary
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]]["BS_data"].extend(promDict[promoter]["BS_data"])
        else:
            #establish promoter count
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]] = {}
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]]["geneInfo"] = promDict[promoter]["geneInfo"]
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]]["geneInfo"]["num_of_promoters"] = 1
            #create "BS_data" dict in geneDict[genename]
            geneDict[promDict[promoter]["geneInfo"]["GeneName"]]["BS_data"] = promDict[promoter]["BS_data"]
    return geneDict

def get_prom_dict_with_QS(threshold = 0):
    # adds QSs to the prom dict in additional category: "QSs"
    # keys are: 0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"
    keyList = [0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum", "BS_sum_all_motifs"]
    promDictQS = get_prom_BS_sum_dict(threshold)
    for prom in promDictQS:
        for key in keyList:
            promDictQS[prom]["QSs"] = {}
            promDictQS[prom]["QSs"][key] = None
        #create Dicts
        sumDict = {"RELA": 0, "RELB": 0, "REL": 0, "all":0}
        avgScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgTSSDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_0 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_1 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_2 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_3 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest15 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest10 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        promDictQS[prom]["avgHighestX"] = {}
        promDictQS[prom]["avgHighestXmsp"] = {}
        avgHighestList = []
        RELA_cW_sum = 0 
        RELA_cW_sur_fit_sum = 0
        for i in range(2,16):
            promDictQS[prom]["avgHighestXmsp"]["avg" + str(i)] = {"RELA": [], "RELB": [], "REL": [], "all": []}
        for i, BS in enumerate(promDictQS[prom]["BS_data"]):
            sumDict["all"] += 1
            sumDict[BS["CRE"]] += 1
            avgScoreDict["all"].append(int(BS["Score"]))
            avgScoreDict[BS["CRE"]].append(int(BS["Score"]))
            avgTSSDict["all"].append(int(BS["dist_TSS"]))
            avgTSSDict[BS["CRE"]].append(int(BS["dist_TSS"]))
            RELA_cW_sum += int(BS["RELA_central_BS_W"]) 
            RELA_cW_sur_fit_sum += int(BS["RELA_cW_sur_fit"])
        # calculate QSs
        for key in avgScoreDict:
            currentScoreList = sorted(avgScoreDict[key], reverse = True)
            if currentScoreList == []:
                QS_0[key] = 0
                QS_1[key] = 0
                QS_2[key] = 0
                QS_3[key] = 0
                avgHighest15[key] = 0
                avgHighest10[key] = 0
                for i in range(2,16):
                    promDictQS[prom]["avgHighestX"]["avg" + str(i)] = 0
                    promDictQS[prom]["avgHighestXmsp"]["avg" + str(i)][key] = 0
            else:
                QS_0[key] = float(round(((mean(currentScoreList) * len(currentScoreList))), 3))
                QS_1[key] = 0
                QS_2[key] = 0
                for rank, score in enumerate(currentScoreList):
                    QS_1[key] += score/(rank + 1)
                    QS_2[key] += score**(1/(rank + 1))
                QS_1[key] = QS_1[key]
                QS_2[key] = QS_2[key]
                QS_3[key] = currentScoreList[0]
                for i in range(2, 16):
                    promDictQS[prom]["avgHighestXmsp"]["avg" + str(i)][key] = mean(currentScoreList[:i])
                if key == "all": # Has to be there so the RelA top BSs dont get used for calculation
                    for i in range(2, 16):
                        avgHighestList.append(mean(currentScoreList[:i]))
                avgHighest15[key] = mean(currentScoreList[:15])
                avgHighest10[key] = mean(currentScoreList[:10])
        if round(mean(avgScoreDict["all"])) == 0.0 or avgScoreDict["all"] == []:
            avgScore= 0
        else:
            avgScore= (mean(avgScoreDict[key]))
        QSList = [QS_0["all"], QS_1["all"], QS_2["all"], QS_3["all"], RELA_cW_sum, RELA_cW_sur_fit_sum, avgScore, sumDict["all"], sumDict]
        for index, QS in enumerate(keyList):
            promDictQS[prom]["QSs"][QS] = QSList[index]
        for i in range(2,16):
            try:
                promDictQS[prom]["avgHighestX"]["avg" + str(i)] = avgHighestList[i-2]
            except IndexError:
                promDictQS[prom]["avgHighestX"]["avg" + str(i)] = 0
                pass
    print("Gene dictionary collected!")
    return promDictQS

def get_gene_dict_with_QS(threshold = 0):
    # adds QSs to the gene dict in additional category: "QSs"
    # keys are: 0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"
    # avgHighestXmsp holds the QS_avg_x_highest_BS for RelA, RelB and c-Rel
    keyList = [0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum","BS_sum_all_motifs"]
    geneDictQS = get_gene_BS_sum_dict(threshold)
    for gene in geneDictQS:
        for key in keyList:
            geneDictQS[gene]["QSs"] = {}
            geneDictQS[gene]["QSs"][key] = None
        #create Dicts
        sumDict = {"RELA": 0, "RELB": 0, "REL": 0, "all":0}
        avgScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgTSSDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_0 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_1 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_2 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_3 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        num_of_promoters = geneDictQS[gene]["geneInfo"]["num_of_promoters"]
        avgHighest15 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest10 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        geneDictQS[gene]["avgHighestX"] = {}
        geneDictQS[gene]["avgHighestXmsp"] = {}
        avgHighestList = []
        RELA_cW_sum = 0 
        RELA_cW_sur_fit_sum = 0
        for i in range(2,16):
            geneDictQS[gene]["avgHighestXmsp"]["avg" + str(i)] = {"RELA": [], "RELB": [], "REL": [], "all": []}
        for i, BS in enumerate(geneDictQS[gene]["BS_data"]):
            sumDict["all"] += 1
            sumDict[BS["CRE"]] += 1
            avgScoreDict["all"].append(int(BS["Score"]))
            avgScoreDict[BS["CRE"]].append(int(BS["Score"]))
            avgTSSDict["all"].append(int(BS["dist_TSS"]))
            avgTSSDict[BS["CRE"]].append(int(BS["dist_TSS"]))
            RELA_cW_sum += int(BS["RELA_central_BS_W"]) 
            RELA_cW_sur_fit_sum += int(BS["RELA_cW_sur_fit"])
        # calculate QSs
        for key in avgScoreDict:
            currentScoreList = sorted(avgScoreDict[key], reverse = True)
            if currentScoreList == []:
                QS_0[key] = 0
                QS_1[key] = 0
                QS_2[key] = 0
                QS_3[key] = 0
                avgHighest15[key] = 0
                avgHighest10[key] = 0
                for i in range(2,16):
                    geneDictQS[gene]["avgHighestX"]["avg" + str(i)] = 0
                    geneDictQS[gene]["avgHighestXmsp"]["avg" + str(i)][key] = 0
            else:
                QS_0[key] = float(round(((mean(currentScoreList) * len(currentScoreList))), 3))
                QS_1[key] = 0
                QS_2[key] = 0
                for rank, score in enumerate(currentScoreList):
                    QS_1[key] += score/(rank + 1)
                    QS_2[key] += score**(1/(rank + 1))
                QS_1[key] = QS_1[key]
                QS_2[key] = QS_2[key]
                QS_3[key] = currentScoreList[0]
                for i in range(2, 16):
                    geneDictQS[gene]["avgHighestXmsp"]["avg" + str(i)][key] = mean(currentScoreList[:i])
                if key == "all": # Has to be there so the RelA top BSs dont get used for calculation
                    for i in range(2, 16):
                        avgHighestList.append(mean(currentScoreList[:i]))
                avgHighest15[key] = mean(currentScoreList[:15])
                avgHighest10[key] = mean(currentScoreList[:10])
        if round(mean(avgScoreDict["all"])) == 0.0 or avgScoreDict["all"] == []:
            avgScore= 0
        else:
            avgScore= (mean(avgScoreDict[key]))
        QSList = [QS_0["all"], QS_1["all"], QS_2["all"], QS_3["all"], RELA_cW_sum, RELA_cW_sur_fit_sum, avgScore, sumDict["all"], sumDict]
        for index, QS in enumerate(keyList):
            geneDictQS[gene]["QSs"][QS] = QSList[index]
        for i in range(2,16):
            try:
                geneDictQS[gene]["avgHighestX"]["avg" + str(i)] = avgHighestList[i-2]
            except IndexError:
                geneDictQS[gene]["avgHighestX"]["avg" + str(i)] = 0
                pass
    print("Gene dictionary collected!")
    return geneDictQS

def get_prom_dict_with_best_QS():
    # the subdictionary with all final QSs as in table 4 (BT), based on the dataset version which led to 
    # the highest odds ratio gets added 
    # holding all the QSs as found in the table: avg4, 1,2,3,0, numBS, avgScore,wang1,wang5,numProm
    # the applied threshold is the same as in the table:  --, 370,310,350,450,480,330,480,--,--
    promDictBestQS = get_prom_dict_with_QS()
    for prom in promDictBestQS:
        promDictBestQS[prom]["best_QS"] = {
                "avg10" : promDictBestQS[prom]["avgHighestX"]["avg10"],
                1 : None,
                2 : None,
                3 : None,
                0 : None,
                "BS_sum" : None,
                "avgScore" : None,
                "Wang1" : None,
                "Wang5" : promDictBestQS[prom]["QSs"]["Wang5"]
        }
    keyList = [1, 2, 3, 0, "BS_sum", "avgScore"]
    for i, threshold in enumerate([370, 310, 350, 450, 480, 330]):
        print(str(i + 1) + " of 7 done :)")
        currentGeneDict = get_prom_dict_with_QS(threshold)
        for prom in promDictBestQS:
            if threshold == 480:
                promDictBestQS[prom]["best_QS"]["Wang1"] = currentGeneDict[prom]["QSs"]["Wang1"]
            promDictBestQS[prom]["best_QS"][keyList[i]] = currentGeneDict[prom]["QSs"][keyList[i]]
    return promDictBestQS

def get_gene_dict_with_best_QS():
    # the subdictionary with all final QSs as in table 4 (BT) gets added 
    # holding all the QSs as found in the table: avg4, 1,2,3,0, numBS, avgScore,wang1,wang5,numProm
    # the applied threshold is the same as in the table:  --, 370,310,350,450,480,330,480,--,--
    geneDictBestQS = get_gene_dict_with_QS()
    for gene in geneDictBestQS:
        geneDictBestQS[gene]["best_QS"] = {
                "avg10" : geneDictBestQS[gene]["avgHighestX"]["avg10"],
                1 : None,
                2 : None,
                3 : None,
                0 : None,
                "BS_sum" : None,
                "avgScore" : None,
                "Wang1" : None,
                "Wang5" : geneDictBestQS[gene]["QSs"]["Wang5"],
                "numProm" : geneDictBestQS[gene]["geneInfo"]["num_of_promoters"]
        }
    keyList = [1, 2, 3, 0, "BS_sum", "avgScore"]
    for i, threshold in enumerate([370, 310, 350, 450, 480, 330]):
        print(str(i + 1) + " of 7 done :)")
        currentGeneDict = get_gene_dict_with_QS(threshold)
        for gene in geneDictBestQS:
            if threshold == 480:
                geneDictBestQS[gene]["best_QS"]["Wang1"] = currentGeneDict[gene]["QSs"]["Wang1"]
            geneDictBestQS[gene]["best_QS"][keyList[i]] = currentGeneDict[gene]["QSs"][keyList[i]]
    return geneDictBestQS

def get_motif_specific_prom_dict(threshold = 0):
    # adds QSs to the prom dict in additional category: "QSs"
    # keys are: 0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"
    keyList = [0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"]
    promDictQS = get_prom_BS_sum_dict(threshold)
    for prom in promDictQS:
        for key in keyList:
            promDictQS[prom]["QSs"] = {}
            promDictQS[prom]["QSs"][key] = None
        #create Dicts
        sumDict = {"RELA": 0, "RELB": 0, "REL": 0, "all":0}
        avgScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgScore = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgTSSDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_0 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_1 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_2 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_3 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest15 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest10 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        promDictQS[prom]["avgHighestX"] = {}
        for i in range(2,16):
            promDictQS[prom]["avgHighestX"]["avg" + str(i)] = {"RELA": 0, "RELB": 0, "REL": 0, "all": 0}
        RELA_cW_sum = 0 
        RELA_cW_sur_fit_sum = 0
        for i, BS in enumerate(promDictQS[prom]["BS_data"]):
            sumDict["all"] += 1
            sumDict[BS["CRE"]] += 1
            avgScoreDict["all"].append(int(BS["Score"]))
            avgScoreDict[BS["CRE"]].append(int(BS["Score"]))
            avgTSSDict["all"].append(int(BS["dist_TSS"]))
            avgTSSDict[BS["CRE"]].append(int(BS["dist_TSS"]))
            RELA_cW_sum += int(BS["RELA_central_BS_W"]) 
            RELA_cW_sur_fit_sum += int(BS["RELA_cW_sur_fit"])
        # calculate QSs
        for key in avgScoreDict:
            currentScoreList = sorted(avgScoreDict[key], reverse = True)
            if currentScoreList == []:
                QS_0[key] = 0
                QS_1[key] = 0
                QS_2[key] = 0
                QS_3[key] = 0
                avgHighest10[key] = 0
            else:
                QS_0[key] = float(round(((mean(currentScoreList) * len(currentScoreList))), 3))
                QS_1[key] = 0
                QS_2[key] = 0
                for rank, score in enumerate(currentScoreList):
                    QS_1[key] += score/(rank + 1)
                    QS_2[key] += score**(1/(rank + 1))
                QS_1[key] = QS_1[key]
                QS_2[key] = QS_2[key]
                QS_3[key] = currentScoreList[0]
                for i in range(2, 16):
                    promDictQS[prom]["avgHighestX"]["avg" + str(i)][key] = mean(currentScoreList[:i])
                avgHighest10[key] = mean(currentScoreList[:10])
            if round(mean(avgScoreDict["all"])) == 0.0 or avgScoreDict["all"] == []:
                avgScore[key] = 0
            else:
                avgScore[key] = (mean(avgScoreDict[key]))
        QSList = [QS_0, QS_1, QS_2, QS_3, RELA_cW_sum, RELA_cW_sur_fit_sum, avgScore, sumDict]
        for index, QS in enumerate(keyList):
            promDictQS[prom]["QSs"][QS] = QSList[index]
    return promDictQS
    
def get_motif_specific_prom_dict_with_best_QS():
    # adds motif specific QSs with optimal dataset versions
    # these can be accessed through the "best_QS" key
    promDictQS = get_motif_specific_prom_dict()
    for prom in promDictQS:
        for key in {"RELA": [], "RELB": [], "REL": [], "all": []}:
            promDictQS[prom]["best_QS"] = {
                "avg10" : promDictQS[prom]["avgHighestX"]["avg10"],
                1 : None,
                2 : None,
                3 : None,
                0 : None,
                "BS_sum" : None,
                "avgScore" : None,
                "Wang1" : None,
                "Wang5" : promDictQS[prom]["QSs"]["Wang5"]
                }
    keyList = [1, 2, 3, 0, "BS_sum", "avgScore"]
    for i, threshold in enumerate([370, 310, 350, 450, 480, 330]):
        print(str(i + 1) + " of 7 done")
        currentPromDict = get_motif_specific_prom_dict(threshold)
        for prom in promDictQS:
            if threshold == 480:
                promDictQS[prom]["best_QS"]["Wang1"] = currentPromDict[prom]["QSs"]["Wang1"]
            promDictQS[prom]["best_QS"][keyList[i]] = currentPromDict[prom]["QSs"][keyList[i]]
    return promDictQS
    
    
def get_motif_specific_gene_dict(threshold = 0):
    # adds QSs to the prom dict in additional category: "QSs"
    # keys are: 0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"
    keyList = [0, 1, 2, 3, "Wang1", "Wang5", "avgScore", "BS_sum"]
    geneDictQS = get_gene_BS_sum_dict(threshold)
    for gene in geneDictQS:
        for key in keyList:
            geneDictQS[gene]["QSs"] = {}
            geneDictQS[gene]["QSs"][key] = None
        #create Dicts
        sumDict = {"RELA": 0, "RELB": 0, "REL": 0, "all":0}
        avgScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgScore = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgTSSDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_0 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_1 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_2 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        QS_3 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest15 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        avgHighest10 = {"RELA": [], "RELB": [], "REL": [], "all": []}
        geneDictQS[gene]["avgHighestX"] = {}
        for i in range(2,16):
            geneDictQS[gene]["avgHighestX"]["avg" + str(i)] = {"RELA": 0, "RELB": 0, "REL": 0, "all": 0}
        RELA_cW_sum = 0 
        RELA_cW_sur_fit_sum = 0
        for i, BS in enumerate(geneDictQS[gene]["BS_data"]):
            sumDict["all"] += 1
            sumDict[BS["CRE"]] += 1
            avgScoreDict["all"].append(int(BS["Score"]))
            avgScoreDict[BS["CRE"]].append(int(BS["Score"]))
            avgTSSDict["all"].append(int(BS["dist_TSS"]))
            avgTSSDict[BS["CRE"]].append(int(BS["dist_TSS"]))
            RELA_cW_sum += int(BS["RELA_central_BS_W"]) 
            RELA_cW_sur_fit_sum += int(BS["RELA_cW_sur_fit"])
        # calculate QSs
        for key in avgScoreDict:
            currentScoreList = sorted(avgScoreDict[key], reverse = True)
            if currentScoreList == []:
                QS_0[key] = 0
                QS_1[key] = 0
                QS_2[key] = 0
                QS_3[key] = 0
                avgHighest10[key] = 0
            else:
                QS_0[key] = float(round(((mean(currentScoreList) * len(currentScoreList))), 3))
                QS_1[key] = 0
                QS_2[key] = 0
                for rank, score in enumerate(currentScoreList):
                    QS_1[key] += score/(rank + 1)
                    QS_2[key] += score**(1/(rank + 1))
                QS_1[key] = QS_1[key]
                QS_2[key] = QS_2[key]
                QS_3[key] = currentScoreList[0]
                for i in range(2, 16):
                    geneDictQS[gene]["avgHighestX"]["avg" + str(i)][key] = mean(currentScoreList[:i])
                avgHighest10[key] = mean(currentScoreList[:10])
            if round(mean(avgScoreDict["all"])) == 0.0 or avgScoreDict["all"] == []:
                avgScore[key] = 0
            else:
                avgScore[key] = (mean(avgScoreDict[key]))
        QSList = [QS_0, QS_1, QS_2, QS_3, RELA_cW_sum, RELA_cW_sur_fit_sum, avgScore, sumDict]
        for index, QS in enumerate(keyList):
            geneDictQS[gene]["QSs"][QS] = QSList[index]
    return geneDictQS

def get_motif_specific_gene_dict_with_best_QS():
    # adds motif specific QSs with optimal dataset versions as foundation
    # these can be accessed through the "best_QS" key
    geneDictQS = get_motif_specific_gene_dict()
    for gene in geneDictQS:
        for key in {"RELA": [], "RELB": [], "REL": [], "all": []}:
            geneDictQS[gene]["best_QS"] = {
                "avg10" : geneDictQS[gene]["avgHighestX"]["avg10"],
                1 : None,
                2 : None,
                3 : None,
                0 : None,
                "BS_sum" : None,
                "avgScore" : None,
                "Wang1" : None,
                "Wang5" : geneDictQS[gene]["QSs"]["Wang5"]
                }
    keyList = [1, 2, 3, 0, "BS_sum", "avgScore"]
    for i, threshold in enumerate([370, 310, 350, 450, 480, 330]):
        print(str(i + 1) + " of 7 done")
        currentGeneDict = get_motif_specific_gene_dict(threshold)
        for gene in geneDictQS:
            if threshold == 480:
                geneDictQS[gene]["best_QS"]["Wang1"] = currentGeneDict[gene]["QSs"]["Wang1"]
            geneDictQS[gene]["best_QS"][keyList[i]] = currentGeneDict[gene]["QSs"][keyList[i]]
    return geneDictQS
    
get_prom_BS_sum_dict()