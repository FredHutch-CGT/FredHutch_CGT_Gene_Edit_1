#!/usr/bin/python3
import sys
import math
import re
import collections
from collections import defaultdict

#
# easy gene editing : read SAM file after genomix alignmend and use CIGAR String 
# to parse and sort results 
#
# Lauren and Mark   3-19-18
# Edit: Lauren 5-14-18
#
cigToken = collections.namedtuple('cigToken',['cmd','position','rest'])
CigEntry = collections.namedtuple('CigEntry', ['cmd', 'start','end','seq'], verbose=False)

#--------------------------------------------------------------------------------
#
# parseInsertCigar - given cigar string cell cigFindNextToken to break it
# into commands of form ('M', start, end)
#
#--------------------------------------------------------------------------------
def cigFindNextToken(cig):
    #
    # look for integers then a single (M,I,D)
    #
    r = re.compile("([0-9]+)([a-zA-Z]+)(.+)")
    m = r.match(cig)
    if m:
        return(cigToken(m.group(2),int(m.group(1)),m.group(3)))
    else:
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(cig)
        if m:
            return(cigToken(m.group(2),int(m.group(1)),""))
        else:
            return(cigToken(False,0,""))
#--------------------------------------------------------------------------------
#
# parseInsertCigar - given cigar string cell cigFindNextToken to break it
# into commands of form ('M', start, end)
#
#--------------------------------------------------------------------------------
def parseInsertCigar(cig):
    #
    # 186M8I67M
    #
    # like to return  (M 0, 254,seq) (I 254 266 NULL) (M 266 300 NULL)
    #
    #(186,M) (8,I) (67,M)
    #Each cigar will have cigEntrys = the number of "commands" in the cigar
    position = 0
    retList = []
    #print("Start ",cig)
    while True:
        ct = cigFindNextToken(cig)
        if (ct.cmd == False):
            return retList
        else:
            if (ct.cmd == "D"):
                newPosition = position
                # keep track of number of deleted base pairs for given command
                numDeletions = ct.position
            else:
                newPosition = position + ct.position
                numDeletions = 0

            cigEntry = CigEntry(ct.cmd,position,newPosition,numDeletions)
            retList.append(cigEntry)
            position = newPosition
            
            cig = ct.rest
            if cig == "":
                return retList

    return False
#--------------------------------------------------------------------------------
#
# easy gene editing : main routine start
#
#--------------------------------------------------------------------------------
f1 = open(sys.argv[1],'r')
#
#  read list of filenames
#
#
# read in sam, only save cigar
#
#
cigIns = {}
#
#
totalSeq = 0
#
#
hdrDict = defaultdict(int)
#
#
for line in f1:
    line = line.strip()
    line = line.split("\t")
    if (line[0][0] != "@"):
        cigar = line[5]
        #
        # deletion only or insertion/deletion
        #
        seq = line[9]
        newCig = parseInsertCigar(cigar)
        #
        # pCig [(cmd,start,end,seq)...()]
        #
        # add the original sequence of insertions to end of cigar string
        # to use as unique key
        #
        parseCig = []
        for entry in newCig:
            if entry.cmd == 'I':
                cigar = cigar + "_" + seq[entry.start:entry.end]
                insSeq = seq[entry.start:entry.end]
                nEntry = CigEntry(entry.cmd,entry.start,entry.end,insSeq)
                parseCig.append(nEntry)
            else:
                parseCig.append(entry)
        #
        # 
        #
        try:
            # cig already in dict
            (tc,tCig) = cigIns[cigar]
            tc += 1
            cigIns[cigar] = (tc,tCig)
        except KeyError:
            #cig is new
            cigIns[cigar] = (1,parseCig)
        #
        # hdr counts
        #
        index = seq.find("GCGGCCGC")
        if index != -1:
          hdrDict[cigar] += 1

#
# sort cigar strings based on count
#
cigSortI = sorted(cigIns.keys(), key=lambda x: cigIns[x][0], reverse=True)
#
# calc number of sequences
#
totalSeq       = 0
totalInsert    = 0
totalDelete    = 0
for (key,tup) in cigIns.items():
    seqCount = tup[0]
    pCig     = tup[1]
    totalSeq += seqCount
    if (key.find('I') != -1):
        totalInsert += 1
    elif (key.find('D') != -1):
        totalDelete += 1
print("----------------------------------------------------------------")
print("Number of sequences                  = ",totalSeq)
print("Number of unique deletion sequences  = ",totalDelete)
print("Number of unique insertion sequences = ",totalInsert)


#
# print out
#
for i in  cigSortI:
    (tc,pCig) = cigIns[i]
    freq = float(tc) / float(totalSeq)
    if (freq > 0.00001):
        print(i,tc,freq)
#
# hdr
#
print("HDR")
allHdr = 0
sort1 =  sorted(hdrDict.items(),key=lambda l: l[1],reverse=True)
for k,v in sort1:
  print(k,v)
  allHdr += v
print("Total HDR = ",allHdr)
#
# write out intermediate file
#
if len(sys.argv) > 3:
    f2 = open(sys.argv[2],'w')
    # last argument should be original amplicon seq for comparison
    orig_seq = sys.argv[3]
    
    f2.write("cigar\tcount\tfreq\tseq\n")
    for i in  cigSortI:
        (tc,pCig) = cigIns[i]
        f2.write(i)
        f2.write("\t")
        f2.write(str(tc))
        f2.write("\t")
        f2.write("{}".format(float(tc)/float(totalSeq)))
        f2.write("\t")
        # counter to keep track of all prior deletions 
        pDelete = 0  
        # loop through pCig commands to rebuild sequences      
        for j in pCig:
            if (j.cmd == "M"):
            	f2.write(orig_seq[j.start+pDelete:j.end+pDelete].lower())
            elif (j.cmd == "D"):
                pDelete += int(j.seq)
                for i in range(1, int(j.seq)+1):
                    f2.write('-')
            elif (j.cmd == "I"):
                f2.write(j.seq.upper())
        f2.write("\n")
    f2.close() 
