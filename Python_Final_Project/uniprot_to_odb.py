#!/usr/bin/env python3

import sys

query_file_name = sys.argv[1]

query_list = list()

ID_to_odb = dict()
odb_to_OGset = dict()
tax_dict = dict()
OG_to_odbList = dict()
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

    for odbID in odb_to_OGset:
        OG_list = odb_to_OGset[odbID]
        OG_set = set(OG_list)
        odb_to_OGset[odbID] = OG_set

print('output of odb_to_OGset:',odb_to_OGset['9606:00232f'])

# the OG_dict is reversed such that a key of an OG will return a list of odb IDs that are in that group
        
for odbID in odb_to_OGlist:
    OG_list = odb_to_OGlist[odbID]

    for ogID in OG_list:
        if ogID in OG_to_odbList:
            OG_to_odbList[ogID].add(odbID)
        else :
            OG_to_odbList[ogID] = set()
            OG_to_odbList[ogID].add(odbID)

print('output of OG_to_odb_list is:',OG_to_odbList['EOG0903023G'])

# a list of OGs associated with an odb will be collected, the list of odbs associated with those OGs are added to a set
# so that a list of odbs is returned for a provided odb

for uniprot_id in query_list:
    homologs = set()
    odb = ID_to_odb[uniprot_id]
    OGlist = odb_to_OGlist[odb]

    for OG in OGlist:
        odbList = OG_to_odbList[OG]

        for related_odb in odbList:
            homologs.add(related_odb)

    uniprot_to_homologs[uniprot_id] = homologs

print(uniprot_to_homologs)
    
# a lookup table is created from the parsed gene dict in which a odb ID as key will return tax ID

#for key in ID_dict:
#    odb = ID_dict[key]['odb']
#    tax = ID_dict[key]['tax']
#    tax_dict[odb] = tax
    
    

#odbID = ID_dict[query]
#OG = OG_dict[odbID]
#tax = tax_dict[OG]

#print(odbID)
#print(OG)
#print(tax)
