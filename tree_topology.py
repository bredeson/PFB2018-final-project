#!/usr/bin/env python

#This script will take a list of taxid and produce a newick tree file of taxids
#and a newick tree file of first common ancestor ancestor(FCA) subtree
#It calculates the number of steps until protein was gained (between FCA & cellular organism root)
#The tree files can be uploaded to itol for visualization and converting taxids to names for you
from ete3 import Tree
from ete3 import NCBITaxa
from ete3 import PhyloTree

ncbi = NCBITaxa()
tree = ncbi.get_topology([7227,9031,9606,562,10090])

tree.write(format=3, outfile="NCBITree_topology.nwk")

#I will need to save it in command line by capturing stdnout, this can be modified to be written in file (optional)
#print(tree.get_ascii(attributes=["sci_name"]))


#get subtree of first common ancestor
first_common_ancestor_subtree = tree.get_common_ancestor("10090", "9031")
first_common_ancestor_subtree.write(format=3, outfile="ancestor_subtree.nwk")
#print(first_common_ancestor_subtree)

#This prints the FCA taxid
first_common_ancestor_taxid = tree.get_common_ancestor("10090", "9031")
print('first_common_ancestor_taxid',first_common_ancestor_taxid.name)

# number of steps till the root (cellular organism) I will need to modify this to be till eukaryote
lineage_list = ncbi.get_lineage(9031)

#This prints the FCA taxon name
first_common_ancestor_name_dict = ncbi.get_taxid_translator(lineage_list)
print(first_common_ancestor_name_dict[int(first_common_ancestor_taxid.name)])
#print(lineage_list)

# -1 is needed to exclude the leaf of interest as it is reportred in the list
# this gives you the number of taxa above the leaf of interest to the cellular organism root
#print(len(lineage_list)-1)

# this gives you the number of taxa above the leaf of interest to the Eukaryota (2759)
print('number of taxa above the leaf of interest to the Eukaryota',len(lineage_list)-3)
# number of steps from the root (cellular organism in this case) to first_common_ancestor
print('number of steps from the root',lineage_list.index(int(first_common_ancestor_taxid.name)))
# number of steps from Eukaryota to first_common_ancestor
print('number of steps from Eukaryota to first_common_ancestor',lineage_list.index(int(first_common_ancestor_taxid.name))-2)
# number of steps from leaf to the FCA
print('number of steps from leaf to the FCA',len(lineage_list)-1-lineage_list.index(int(first_common_ancestor_taxid.name)))


#Doing the same for the other leaf of interest
lineage_list = ncbi.get_lineage(10090)
# this gives you the number of taxa above the leaf of interest to the cellular organism root
print('number of taxa above the leaf of interest',len(lineage_list)-1)

# number of steps from other leaf to the FCA
#print(len(lineage_list)-1-lineage_list.index(int(first_common_ancestor_taxid.name)))
