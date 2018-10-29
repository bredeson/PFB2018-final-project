#!/usr/bin/env python

#This script will take tab-separated uniprotid and taxids and finds first_recenet_common ancestor
#and number of steps from FCA to Eukaryota

#This script uses ete3 for dealing with trees and graphs
import sys
from ete3 import Tree
from ete3 import NCBITaxa
from ete3 import PhyloTree
ncbi = NCBITaxa()
#Here I open the file that Matt script creates and loops in each line and get the taxids
with open('query_list_output2.tab','r') as fo:
    for line in fo:
        line = line.rstrip()
        (uniprotid,taxids) = line.split('\t')
        one_taxid = taxids.split(',') # divide the list of taxids to diff taxids 'strings'
        for i in range(len(one_taxid)):
            one_taxid[i] = int(one_taxid[i]) #ete3 take a list of taxid integers
        #print(one_taxid)
        tree = ncbi.get_topology(one_taxid) #creates the tree of taxids for each uniprot id
        for i in range(len(one_taxid)):
            one_taxid[i] = str(one_taxid[i]) # returning to list of strings again
        #print(one_taxid)
        for i in range(len(one_taxid)):
            one_taxid_str = one_taxid[i] # creating separate strings from taxids list of strings
        #get the first common ancestor of each uniprotid as taxid
        first_common_ancestor_taxid = tree.get_common_ancestor(one_taxid_str)
        #print('first_common_ancestor_taxid',first_common_ancestor_taxid.name)

        #get a lineage of one taxid
        lineage_list = ncbi.get_lineage(one_taxid_str)
        #This prints the FCA taxon name
        #translate the lineage to names and get FCA as a taxon name
        first_common_ancestor_name_dict = ncbi.get_taxid_translator(lineage_list)
        #print(first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)])
        #print('number of steps from Eukaryota to first_common_ancestor',lineage_list.index(int(first_common_ancestor_taxid.name))-2)

        #prints all the results
        print(uniprotid,first_common_ancestor_taxid.name,first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)],
        lineage_list.index(int(first_common_ancestor_taxid.name))-2,tree.write(format=3))
