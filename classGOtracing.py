#!/usr/bin/env python3
import re
import pronto
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from wordcloud import WordCloud
#ont = pronto.Ontology("/projects/jessens-angels/sla/PFB2018-final-project/go.owl")

class GOclass:
    def __init__ (self, owlFile):
        self.ont = pronto.Ontology(owlFile)


    def uniprot2GO(self, uniprotID, gafFile): #gets nonredundant set of uniprot IDs for given uniprot ID as a string
        GO_id = set()
        GO_annot = []
        with open(gafFile, "r") as file_object:
            for line in file_object:
                if uniprotID in line:
                   found =  re.findall(uniprotID+'\s[\w\-\s]+(GO\:\d+)\s', line)
                   if found != []:
                       GO_id.add(found[0])
        for goID in GO_id:
            GO_annot.append(self.ont[goID].name)

        return GO_annot 

    def freqGOannot(self, GO_annot): #submit list of GO annotations
        GOfreqDict = {}
        for term in GO_annot:
            GOfreqDict[term] = GO_annot.count(term)
        return GOfreqDict

    def plotGOannot(self, GofreqDict, filename): #takes a dictionary of GOterms and frequencies, and a filename to write in, and generates wordcloud of GO terms, writes in file
        wordcloud = WordCloud()
        wordcloud.fit_words(GofreqDict)
        plt.imshow(wordcloud, interpolation = 'bilinear')
        plt.axis("off")
        plt.savefig(filename)

#annot, GOlist = uniprot2GO('P26640')

if __name__ == '__main__':
    testOnt = GOclass('/projects/jessens-angels/sla/PFB2018-final-project/go.owl')
    print(testOnt.uniprot2GO('P26640', "goa_human.gaf"))

#plotGOannot(dictionary,"testfig.png")
