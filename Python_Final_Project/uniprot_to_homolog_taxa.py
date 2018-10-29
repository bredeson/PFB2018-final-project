#!/usr/bin/env python3

import sys

query_file_name = sys.argv[1]
output_file_name = sys.argv[2]

query_list = list()

ID_to_odb = dict()
odb_to_OGset = dict()
tax_dict = dict()
OG_to_odbSet = dict()
tax_dict = dict()
uniprot_to_homologs = dict()

# the query file is read in and a list of UniProt IDs is generated

with open(query_file_name, 'r') as file_obj:
    for line in file_obj:
        line = line.rstrip()
        query_list.append(line)

# the genes.tab file is parsed into a two-level dictionary where the odb ID and tax ID are saved as the value for a uniprot ID key
            
with open('odb9v1_genes.tab', 'r') as file_obj:
   for line in file_obj:
        line_list = line.split()
        ID_to_odb[line_list[3]] = {'odb' : line_list[0],'tax': line_list[1]}

print('output of ID_to_odb:',ID_to_odb['O14802']) 

# the OG2genes.tab file is parse into a dictionary where a key of odb ID will return the OG
        
with open('odb9v1_OG2genes.tab', 'r') as file_obj:
    for line in file_obj:
        odbID = line.split()[1]
        OG = line.split()[0]
        if odbID in odb_to_OGset:
            odb_to_OGset[odbID].append(OG)
        else :
            odb_to_OGset[odbID] = [OG]

# this converts the list of OGs associated with each odb to a set
            
    for odbID in odb_to_OGset:
        OG_list = odb_to_OGset[odbID]
        OG_set = set(OG_list)
        odb_to_OGset[odbID] = OG_set

print('output of odb_to_OGset:',odb_to_OGset['9606:00232f'])

# the OG_dict is reversed such that a key of an OG will return a list of odb IDs that are in that group
        
for odbID in odb_to_OGset:
    OG_list = odb_to_OGset[odbID]

    for ogID in OG_list:
        if ogID in OG_to_odbSet:
            OG_to_odbSet[ogID].add(odbID)
        else :
            OG_to_odbSet[ogID] = set()
            OG_to_odbSet[ogID].add(odbID)

print('output of OG_to_odb_list is:',OG_to_odbSet['EOG0903023G'])

# a list of OGs associated with an odb will be collected, the list of odbs associated with those OGs are added to a set
# so that a list of odbs is returned for a provided odb

for uniprot_id in query_list:
    homologs = set()
    try:
        odb = ID_to_odb[uniprot_id]['odb']
    except KeyError:
        print('Can not find',uniprot_id,'in ortho db')
        continue
    try:
        OGlist = odb_to_OGset[odb]
    except KeyError:
        print('Can not find an OG for',odb)
        continue
    for OG in OGlist:
        odbList = OG_to_odbSet[OG]

        for related_odb in odbList:
            homologs.add(related_odb)

    uniprot_to_homologs[uniprot_id] = homologs

# the homologs for each entry in the input file are found and the taxa from these homologs are written to an output file

masterset = set()

with open(output_file_name , 'w') as outputFile:
    for uniprot_ID in uniprot_to_homologs:
        mapped_odbs = uniprot_to_homologs[uniprot_ID]

        mapped_odb_list = []
        
        for odb in mapped_odbs:
            taxa = odb.split(':')[0]
            mapped_odb_list.append(taxa)

        mapped_odb_set = set(mapped_odb_list)
        masterset = masterset.union(mapped_odb_set)
        
        mapped_odb_str = ''

        for item in mapped_odb_set:
            mapped_odb_str += item + ','

        mapped_odb_str = mapped_odb_str.rstrip(',')
        
        outputFile.write(uniprot_ID+'\t'+mapped_odb_str+'\n')
            
with open('all_taxa_output_SP.txt','w') as outMaster:
    masterstr = ''
    for taxa in masterset:
        masterstr += taxa + ','
    masterstr = masterstr.rstrip(',')
    outMaster.write(masterstr)
