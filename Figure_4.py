#creates list of all ensembls present in both bielefeld as well as the other datasets
#filename is "bi_all_overlap_ensembls.csv"
#prints the data for figure 6
#ATTENTION: This script requires the geckodriver binary to be in in a folder named geckodriver in the scripts directory
# a manual step was done: The genenames were copied from the Yang et al., 2016 paper and were updated with the ncbi symbol checker (link in workflow) and saved as "symbol_checker_results_seu.txt" in the data folder. Same was done for the Bielefeld dataset(see line 274f)

# Workflow:
# Boston and Lille 
# -> extract and update RefSeqs using selenium, 
# -> use mapping File to get Ensembl IDs for all NF-kB target Genes
#SEU and Bielefeld database :
# -> update genesymbols found in the paper using: https://www.genenames.org/cgi-bin/symbol_checker
# -> use Rest API and EPD database to add Ensembl IDs to the gene names
# Final:
# Calculate Overlap between the databases and print it
import os
import re
import csv
import sys
import requests
from selenium import webdriver
import time
import lxml.html as html
from selenium.common.exceptions import TimeoutException
import platform

def openReadableCSV(filename):
    readableCSV=csv.reader(open(str(filename),newline=''),delimiter='\t')
    return readableCSV

def returnReadableCSVList(filename):
    readableCSV=csv.reader(open(str(filename)),delimiter='\t')
    readableCSVList=list(readableCSV)
    return readableCSVList

def openWriteableCSV(filename):
    writeableCSV=csv.writer(open(str(filename),'w+',newline=''),delimiter='\t')
    return writeableCSV
    
def extractEPDGeneName(PromoterID):
    EPDName=re.compile(r'([a-zA-Z0-9@\-.]+)_\d+').findall(PromoterID)
    if EPDName==[]:
        print("No EPD Name found in "+PromoterID)
        return "Not found"
    else:
        return EPDName[0]
    
# create needed folders if they are not present
if not os.path.exists(os.path.join(".", "data")):
    os.mkdir(os.path.join(".", "data"))
docuFile = openWriteableCSV(os.path.join(".","data", "docuFileFigure6.txt"))
if not os.path.exists(os.path.join(".", "source_files")):
    os.mkdir(os.path.join(".", "source_files"))
    
#check whether lille and boston refs are there from another try 
elif(os.path.exists(os.path.join(".", "source_files", "updatedBostonRefs.txt")) and
    os.path.exists(os.path.join(".", "source_files", "updatedLilleRefs.txt"))):
    bostonRefFile = openReadableCSV(os.path.join(".", "source_files", "updatedBostonRefs.txt"))
    bostonRefsList = []
    for refseq in bostonRefFile:
        bostonRefsList.append(refseq[0])
    lilleRefFile = openReadableCSV(os.path.join(".", "source_files", "updatedLilleRefs.txt"))
    lilleRefsList = []
    for refseq in lilleRefFile:
        lilleRefsList.append(refseq[0])
else:
    #get boston refseqs
    boston = requests.get('https://www.bu.edu/nf-kb/gene-resources/target-genes/')
    if(boston.status_code==200):
        print('Status: Ok!')
    else:
        print('Status: I failed to download the website!')
        input()
    bostonRefsList=re.compile(r'([A-Z]+_*[A-Z]*\d{5,200})\.*\d*').findall(boston.text[18000:])
    print(str(len(bostonRefsList)), " RefSeqs were extracted from the boston database")
    bostonRefsSet = set(bostonRefsList)
    print(str(len(bostonRefsSet)), " unique RefSeqs remain after removing duplicates")
    #get lille refseqs
    lille = requests.get('http://bioinfo.lifl.fr/NF-KB/')
    if(lille.status_code==200):
        print('Status: Ok!')
    else:
        print('Status: I failed to download the website!')
        input()
    lilleRefsList=re.compile(r'target="RefSeq">(.+)</a>').findall(lille.text[3000:100000])
    print(str(len(lilleRefsList)), " RefSeqs were extracted from the Lille database")
    lilleRefsSet = set(lilleRefsList)
    print(str(len(lilleRefsSet)), " unique RefSeqs remain after removing duplicates")
    #update boston and lille refseqs
    def updateRefs(refSeqList, amountOfRefSeqs):
        if platform.system()=="Windows":
            browser = webdriver.Firefox(executable_path=os.path.join('.', 'geckodriver', 'geckodriver.exe'))
        if platform.system()=="Darwin":
            browser = webdriver.Firefox(executable_path=os.path.join('.','geckodriver', 'geckodriver'))
        updatedRefs = 0
        for i in range(amountOfRefSeqs):
            try:
                checkText=browser.get('https://www.ncbi.nlm.nih.gov/nuccore/'+str(refSeqList[i]))
                time.sleep(1)
            except TimeoutException:
                print("timeout")
                time.sleep(3)
            tmpWebsiteContent=html.document_fromstring(browser.page_source).text_content()
            humanityCheck = re.compile(r'SOURCE\s+Homo\ssapiens\s\(human\)').search(tmpWebsiteContent)
            if humanityCheck == None:
                refSeqList[i] = None
                continue
            status=re.compile(r'NCBI\sReference\sSequence:\s('+refSeqList[i]+'.\d+)').search(tmpWebsiteContent)
            if status != None:
                refSeqList[i] = status.group(1)
            else:
                status=re.compile(r'This\ssequence\shas\sbeen\sreplaced\sby\s([A-Z_]+\d+.*\d*)').search(tmpWebsiteContent)
                if(status != None):
                    updatedRefs+=1
                    while(status!=None):
                        refSeqList[i]=status.group(1)
                        try:
                            checkText=browser.get('https://www.ncbi.nlm.nih.gov/nuccore/'+str(refSeqList[i]))
                            time.sleep(1)
                        except TimeoutException:
                            time.sleep(3)
                        tmpWebsiteContent=html.document_fromstring(browser.page_source).text_content()
                        # does not recognize all sequences...
                        status=re.compile(r'This\ssequence\shas\sbeen\sreplaced\sby\s([A-Z]+_\d+.*\d*)').search(tmpWebsiteContent)
                else:
                    status=re.compile(r'No\sitems\sfound.').search(tmpWebsiteContent)
                    if(status!=None):
                        refSeqList[i]=None
                status=re.compile(r'Record\sremoved').search(tmpWebsiteContent)
                if(status!=None):
                    refSeqList[i]=None
        refSeqList = list(filter(lambda x: x != None, refSeqList))
        print(str(updatedRefs), " RefSeqs were updated")
        return refSeqList

    bostonRefsList = updateRefs(list(bostonRefsSet), len(bostonRefsSet))
    print(bostonRefsList)
    lilleRefsList = updateRefs(list(lilleRefsSet), len(lilleRefsSet))
    print(lilleRefsList)
    #write Refseqs to file so the programm runs shorter during testing
    bostonRefFile = openWriteableCSV(os.path.join(".", "source_files", "updatedBostonRefs.txt"))
    for refseq in bostonRefsList:
        bostonRefFile.writerow([refseq])
    lilleRefFile = openWriteableCSV(os.path.join(".", "source_files", "updatedLilleRefs.txt"))
    for refseq in lilleRefsList:
        lilleRefFile.writerow([refseq])
#create refseq to Ensembl dictionary using a mapping file from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
ref2EnsemblFile = openReadableCSV(os.path.join(".", "source_files", "gene2ensembl.txt"))
ref2EnsDict = {}
category = {}
for header in ref2EnsemblFile:
    for index, columnName in enumerate(header):
        category[columnName] = index
    break
for row in ref2EnsemblFile:
    ref2EnsDict[row[category['RNA_nucleotide_accession.version']]] = row[category['Ensembl_gene_identifier']]
#make Lists with boston and lille Ensembls
print("boston not found in ensembl Dict genes: ")
bostonEnsembls = []
for i in bostonRefsList:
    try:
        bostonEnsembls.append(ref2EnsDict[i])
    except:
        print(i)
        continue
print("lille not found in ensembl Dict genes: ")
lilleEnsembls = []
for i in lilleRefsList:
    try:
        lilleEnsembls.append(ref2EnsDict[i])
    except:
        print(i)
        continue
#write Ensembls to file
bostonEnsFile = openWriteableCSV(os.path.join(".", "data", "bostonEnsembls.txt"))
for i in bostonEnsembls:
    bostonEnsFile.writerow([i])
lilleEnsFile = openWriteableCSV(os.path.join(".", "data", "lilleEnsembls.txt"))
for i in lilleEnsembls:
    lilleEnsFile.writerow([i])
#make lists with SEU Ensembls
def getEnsemblFromGeneSymbol(GeneSymbol):
    server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/"+GeneSymbol+"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    if decoded != []:
        restAPIEnsemblSet = set()
        for i in range(len(decoded)):
            restAPIEnsemblSet.add(decoded[0]["id"])
    else:
        restAPIEnsemblSet = set()
    crossRefFile = openReadableCSV(os.path.join(".", "source_files", "cross_references.txt"))
    for row in crossRefFile:
        EPDEnsemblID = None
        if row[1] == GeneSymbol:
            EPDEnsemblID = row[0]
    if EPDEnsemblID in restAPIEnsemblSet:
        return EPDEnsemblID
    elif len(restAPIEnsemblSet) == 1:
        return list(restAPIEnsemblSet)[0]
    elif len(restAPIEnsemblSet) > 1:
        print(restAPIEnsemblSet)
        input()
    else:
        return []
        print("no Ensembl found for", GeneSymbol)
    #get Ensembls for Southeast University Gene
if( os.path.exists(os.path.join(".", "data", "seuEnsembls.txt"))):
    seuEnsemblList = []
    seuEnsemblFile = openReadableCSV(os.path.join(".", "data", "seuEnsembls.txt"))
    for row in seuEnsemblFile:
        seuEnsemblList.append(row[0])
else:
    seu_gNs = returnReadableCSVList(os.path.join(".", "data", "symbol_checker_results_seu.txt"))[1:]
    seuEnsemblFile = openWriteableCSV(os.path.join(".", "data", "seuEnsembls.txt"))
    print("SEU genes extracted from paper: ", str(len(seu_gNs)))
    seuEnsemblList = []
    for index,i in enumerate(seu_gNs):
        print(index)
        if i[1] != "Unmatched":
            Ensembls = ""
            count = 0
            currentEnsemblID = getEnsemblFromGeneSymbol(i[2])
            if currentEnsemblID == []:
                currentEnsemblID = getEnsemblFromGeneSymbol(i[0])
            for j in [currentEnsemblID]:
                if (count+1) == len([currentEnsemblID]) or len([currentEnsemblID]) == 1:
                    Ensembls += j
                else:
                    Ensembls += j + "_%_"
                count += 1
            if count > 1:
                print("more than one!")
                input()
                seuEnsemblFile.writerow(Ensembls.split("_%_"))
                seuEnsemblList.append(Ensembls.split("_%_"))
            elif count == 1:
                seuEnsemblFile.writerow([Ensembls])
                seuEnsemblList.append(Ensembls)
            else:
                print("no Ensembl ID found")
seuEnsemblList = list(filter(lambda x: x != None and x != [] , seuEnsemblList))
seuEnsemblSet = set(filter(lambda x: x != None and x != [] , seuEnsemblList))
print("Amount of SEU genes for which Ensembl IDs where found: ", str(len(seuEnsemblList)))
print("Amount of duplicates: ", str(len(seuEnsemblList) - len(set(seuEnsemblSet))))
seuEnsemblList = list(seuEnsemblSet)

#make list with bielefeld Ensembls
if( os.path.exists(os.path.join(".", "data", "biEnsembls.txt"))):
    biEnsemblList = []
    biEnsemblFile = openReadableCSV(os.path.join(".", "data", "biEnsembls.txt"))
    for row in biEnsemblFile:
        biEnsemblList.append(row[0])
else:
    def extractEPDGeneName(PromoterID):
        EPDName=re.compile(r'([a-zA-Z0-9@\-.]+)_\d+').findall(PromoterID)
        if EPDName==[]:
            print("No EPD Name found in "+PromoterID)
            return "Not found"
        else:
            return EPDName[0]
    biAllData = openReadableCSV(os.path.join(".", "source_files", "all_data.txt"))
    biGeneNameFile = openWriteableCSV(os.path.join(".", "data", "bi_gene_names.txt"))
    biGeneNameSet = set()
    for row in biAllData:
        if biAllData.line_num > 1:
            biGeneNameSet.add(extractEPDGeneName(row[0]))
    for geneName in biGeneNameSet:
        biGeneNameFile.writerow([geneName])
#biGeneNameFile needs to be checked with: https://www.genenames.org/cgi-bin/symbol_checker
#this has to be done manually. The file needs to be named and placed as seen in the next line
    try:
        bi_gNs = returnReadableCSVList(os.path.join(".", "data", "symbol_checker_results_bi.txt"))[1:]
    except:
        print("please use https://www.genenames.org/cgi-bin/symbol_checker on \nthe file 'bi_gene_names.txt' file in the 'data' folder and save the results as\n 'data/symbol_checker_results_bi.txt' in the data folder.")
        print("if you press Enter so an the program will continue")
    finally:
        bi_gNs = returnReadableCSVList(os.path.join(".", "data", "symbol_checker_results_bi.txt"))[1:]

    biEnsemblFile = openWriteableCSV(os.path.join(".", "data", "biEnsembls.txt"))
    print("Number of Bielefeld genes extracted: ", str(len(bi_gNs)))
    biEnsemblList = []
    usedBiGeneNameSet = set()
    typeCount = [0, 0]
    for index,i in enumerate(bi_gNs):
        if index % 1000 == 0:
            print(index)
        if i[1] != "Unmatched" and i[0] not in usedBiGeneNameSet:
            usedBiGeneNameSet.add(i[0])
            Ensembls = ""
            count = 0
            currentEnsemblID = getEnsemblFromGeneSymbol(i[2])
            if currentEnsemblID == []:
                currentEnsemblID = getEnsemblFromGeneSymbol(i[0])
            if currentEnsemblID != []:
                for j in [currentEnsemblID]:
                    if (count+1) == len([currentEnsemblID]) or len([currentEnsemblID]) == 1:
                        Ensembls += j 
                    else:
                        Ensembls += j + "_%_"
                    count += 1
            if count > 1:
                biEnsemblFile.writerow(Ensembls.split("_%_"))
                biEnsemblList.append(Ensembls)
                typeCount[0] += 1
            elif count == 1:
                biEnsemblFile.writerow([Ensembls])
                biEnsemblList.append(Ensembls)
                typeCount[1] += 1
            else:
                print("no Ensembl ID found")
    del usedBiGeneNameSet
    print("types of annotation used: " + str(typeCount))
biEnsemblList = list(filter(lambda x: x != None and x != [] , biEnsemblList))
biEnsemblSet = set(biEnsemblList)
print("Amount of Bielefeld genes for which Ensembl IDs where found: ", str(len(biEnsemblList)))
print("Amount of duplicates: ", str(len(biEnsemblList) - len(set(biEnsemblSet))))
biEnsemblList = list(biEnsemblSet)
print("bi ensembls ready")

#calculate Overlap
def getOverlap(list1, list2, usdEns):
    overlapList=[]
    for i in list1:
        if isinstance(i, list):
            for j in i:
                if j not in usdEns:
                    for k in list2:
                        if isinstance(k, list):
                            if j in k:
                                overlapList.append(j)
                                break
                        elif isinstance(k, str):
                            if j == k:
                                overlapList.append(j)
                                break
                        else:
                            print("check your code! 1")
                            input()             
        elif isinstance(i, str):
            if i not in usdEns:
                for L in list2:
                    if isinstance(L, list):
                        if i in L:
                            overlapList.append(i)
                            break
                    elif isinstance(L, str):
                        if i == L:
                            overlapList.append(i)
                            break
                    else:
                        print("check your code! 2")
                        input()     
        else:
            print("checkyourcode! 3")
            input()
    return overlapList

def getOverlapDict():
    biEnsembls = biEnsemblList
    liEnsembls = lilleEnsembls
    boEnsembls = bostonEnsembls
    seuEnsembls = seuEnsemblList
                    
    ensemblDict = {}
    ensemblDict["li"] = liEnsembls
    print("liEns: "+str(len(liEnsembls)))
    ensemblDict["bo"] = boEnsembls
    print("boEns: "+str(len(boEnsembls)))
    ensemblDict["bi"] = biEnsembls
    print("biEns: "+str(len(biEnsembls)))
    ensemblDict["seu"] = seuEnsembls
    print("seuEns: "+str(len(seuEnsembls)))
    ensemblDict["biLiBoSeu"] = getOverlap(getOverlap(getOverlap(ensemblDict["li"], ensemblDict["bo"], []),
    ensemblDict["bi"], []),
    ensemblDict["seu"], [])
    #print(len(ensemblDict["biLiBoSeu"]))
    global usdEns 
    usdEns = set() # Ensembls which are already used and should be excluded
    usdEns.update(ensemblDict["biLiBoSeu"])

    ensemblDict["liBiBo"] = getOverlap(getOverlap(ensemblDict["bo"], ensemblDict["li"], usdEns),
    ensemblDict["bi"], usdEns)
    usdEns.update(ensemblDict["liBiBo"])

    ensemblDict["liSeuBo"] = getOverlap(getOverlap(ensemblDict["li"], ensemblDict["seu"], usdEns),
    ensemblDict["bo"], usdEns)
    usdEns.update(ensemblDict["liSeuBo"])

    ensemblDict["biBoSeu"] = getOverlap(getOverlap(ensemblDict["bi"], ensemblDict["seu"],usdEns),
    ensemblDict["bo"], usdEns)
    usdEns.update(ensemblDict["biBoSeu"])

    ensemblDict["biLiSeu"] = getOverlap(getOverlap(ensemblDict["li"], ensemblDict["seu"],usdEns),
    ensemblDict["bi"], usdEns)
    usdEns.update(ensemblDict["biLiSeu"])

    ensemblDict["liSeu"] = getOverlap(ensemblDict["li"], ensemblDict["seu"], usdEns)
    usdEns.update(ensemblDict["liSeu"])

    ensemblDict["boSeu"] = getOverlap(ensemblDict["bo"], ensemblDict["seu"], usdEns)
    usdEns.update(ensemblDict["boSeu"])

    ensemblDict["biSeu"] = getOverlap(ensemblDict["bi"], ensemblDict["seu"], usdEns)
    usdEns.update(ensemblDict["biSeu"])

    ensemblDict["liBo"] = getOverlap(ensemblDict["li"], ensemblDict["bo"], usdEns)
    usdEns.update(ensemblDict["liBo"])

    ensemblDict["liBi"] = getOverlap(ensemblDict["li"], ensemblDict["bi"], usdEns)
    usdEns.update(ensemblDict["liBi"])

    ensemblDict["boBi"] = getOverlap(ensemblDict["bo"], ensemblDict["bi"], usdEns)
    usdEns.update(ensemblDict["boBi"])

    return ensemblDict
    
def printVennData():
    ensemblDict = getOverlapDict()
    for i in ["li", "bi", "bo", "seu"]:
        print(i + ": " + str(len(ensemblDict[i])))
        tmpList=[]
        for j in ensemblDict[i]:
            if isinstance(j, list):
                for k in j:
                    if k not in usdEns:
                        tmpList.append(j)
                        break
            elif isinstance(j, str):
                if j not in usdEns:
                    tmpList.append(j)
            else:
                print("problem!")
                input()
        ensemblDict[i] = tmpList
        print(i+str(len(ensemblDict[i])))
    for i in ensemblDict.keys():
        print(str(i) + ": " + str(len(ensemblDict[i])))
    # save ensembls overlapping between bielefeld and other datasets in a file
    with open(os.path.join(".", "source_files", "bi_all_overlap_ensembls.csv"), "w+", newline = "") as f:
        overlapEnsmblFile = csv.writer(f, delimiter = "\t")
        overlapEnsembls = []
        for i in ["biLiBoSeu", "biBoSeu", "biLiSeu", "liBiBo", "biSeu", "liBi", "boBi"]:
            overlapEnsembls.extend(i)
        for i in overlapEnsembls:
            overlapEnsmblFile.writerow([i])

printVennData()    