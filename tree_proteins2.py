#!/usr/bin/env python

#This script will take tab-separated uniprotid and taxids and finds first_recenet_common ancestor
#and number of steps from FCA to Eukaryota

#This script uses ete3 for dealing with trees and graphs
import sys
from ete3 import Tree
from ete3 import NCBITaxa
from ete3 import PhyloTree
ncbi = NCBITaxa()
tax_dict = dict()
tree_dict = dict()
ontology = dict()
#Here I open the file that Matt script creates and loops in each line and get the taxids
with open('query_list_output2.tab','r') as fo:
    for line in fo:
        line = line.rstrip()
        (uniprotid,taxids) = line.split('\t')
        one_taxid = taxids.split(',') # divide the list of taxids to diff taxids 'strings'
        tax_dict[uniprotid] = one_taxid
        one_taxid_int = []
        for i in range(len(one_taxid)):
            one_taxid_int.append(int(one_taxid[i])) #ete3 take a list of taxid integers
        #print(one_taxid)
        tree = ncbi.get_topology(one_taxid) #creates the tree of taxids for each uniprot id
        tree_dict[uniprotid] = tree
        one_taxid_str = []
        for i in range(len(one_taxid_int)):
            one_taxid_str.append(str(one_taxid_int[i])) # returning to list of strings again
        #print(one_taxid)
        for i in range(len(one_taxid_str)):
            one_taxon_id = one_taxid_str[i] # creating separate strings from taxids list of strings
        #get the first common ancestor of each uniprotid as taxid
        first_common_ancestor_taxid = tree.get_tree_root()
        #print(first_common_ancestor_taxid.name)
        if first_common_ancestor_taxid.name in ontology:
            uniprotid_set = ontology[first_common_ancestor_taxid.name]
            uniprotid_set.add(uniprotid)
        else:
            ontology[first_common_ancestor_taxid.name] = set()
            ontology[first_common_ancestor_taxid.name].add(uniprotid)


        #print('first_common_ancestor_taxid',first_common_ancestor_taxid.name)
        #print(first_common_ancestor_taxid)
        #get a lineage of one taxid
        #ontology[first_common_ancestor_taxid] = uniprotid

        lineage_list = ncbi.get_lineage(one_taxon_id)
        #This prints the FCA taxon name
        #translate the lineage to names and get FCA as a taxon name
        first_common_ancestor_name_dict = ncbi.get_taxid_translator(lineage_list)
        #print(first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)])
        #print('number of steps from Eukaryota to first_common_ancestor',lineage_list.index(int(first_common_ancestor_taxid.name))-2)

        #prints all the summary results
        print(uniprotid,first_common_ancestor_taxid.name,first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)],
        lineage_list.index(int(first_common_ancestor_taxid.name))-2,tree.write(format=3))
'''        with open('trees.tab' , 'w') as outputFile:
            outputFile.write(uniprotid+'\t'+tree.write(format=3))'''
        #print(list(tree))
#with open('summary.tab' , 'w') as outputFile:
#    outputFile.write(uniprotid,first_common_ancestor_taxid.name,first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)],
#    lineage_list.index(int(first_common_ancestor_taxid.name))-2,tree.write(format=3))
print(ontology)
with open('all_taxa_output.txt','r') as fo:
    for line in fo:
        line = line.rstrip()
        one_taxid = line.split(',') # divide the list of taxids to diff taxids 'strings'
        #print(one_taxid)
        all_tree = ncbi.get_topology(one_taxid) #creates the tree of taxids for each uniprot id
all_tree.write(format=3, outfile="all_proteins_treezz.nwk")
