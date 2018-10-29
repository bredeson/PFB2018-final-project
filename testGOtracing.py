#!/usr/bin/env python3
from classGOtracing import*
import re
import sys

#sys.argv[1] = owl file, sys.argv[2] = input uniprot file, sys.argv[3] human GOannot GAF file
#generate ontology class
ont = GOclass(sys.argv[1])

#parse uniprot data file
with open(sys.argv[2], "r") as file_object:
    uniprot_dicts = {}
    for line in file_object:
        found = re.findall(r'(\d+)\t\{([\w\'\,\s]+)\}',line)[0] #returns list of tuples, only take the tuple with (taxonid, stringofUniprots)
        uniprot_dicts[found[0]] = found[1] #creates dictionary of {taxonid of last common ancestor: Uniprotstring}

#generate a dictionary, key = taxonid, value = list with all the GO annotations
GO_annot = {}
for ancestor in uniprot_dicts:
    GO_annot[ancestor] = []
    uniprots_cln = uniprot_dicts[ancestor].replace("'", "")
    uniprots_cln2 = uniprots_cln.replace(" ","")
    uniprots = uniprots_cln2.split(",")
    for uniprotID in uniprots:
        GOlist = ont.uniprot2GO(uniprotID,sys.argv[3])
        GO_annot[ancestor] += GOlist
#    print (len(GO_annot[ancestor]))
#print (GO_annot)
#now generate new dictionary, key = taxonid, value = dictionary of frequencies of GOannotations
BgdGOFreq = {} #just a dictionary of frequencies of all terms
GOannotFreq = {} #dictionary of dictionary, frequencies of terms(values) in each taxa (keys)
GOFreq = open("GOfrequencies.txt","w")
bgd = open("GObackgroundfreq.txt","w")
for ancestor in GO_annot:
#    GOannotFreq[ancestor] = {}
    freq = ont.freqGOannot(GO_annot[ancestor])
    GOannotFreq[ancestor] = freq
    for term in GOannotFreq[ancestor]:
        if term in BgdGOFreq:
            BgdGOFreq[term] += GOannotFreq[ancestor][term]
        else:
            BgdGOFreq[term] = GOannotFreq[ancestor][term]
    GOFreq.write('{}\t{}\n'.format(ancestor,GOannotFreq[ancestor]))
for term in BgdGOFreq:
    bgd.write('{}\t{}\n'.format(term, BgdGOFreq[term]))
#for each age - generate a wordcloud named after the taxonID of the ancestor
for ancestor in GOannotFreq:
    filename = ancestor + '.png'
    ont.plotGOannot(GOannotFreq[ancestor], filename)
