#!/usr/bin/env python3

# the query file is read in and a list of UniProt IDs is generated

def query_list_parser(query_file_name):

    query_list = []

    with open(query_file_name, 'r') as file_obj:
        for line in file_obj:
            line = line.rstrip()
            query_list.append(line)

    return query_list

# the genes.tab file is parsed into a two-level dictionary where the odb ID and tax ID are saved as the value for a uniprot ID key

def make_ID_to_odb():

    ID_to_odb = dict()   

    with open('odb9v1_genes.tab', 'r') as file_obj:
        for line in file_obj:
            line_list = line.split()
            ID_to_odb[line_list[3]] = {'odb' : line_list[0],'tax': line_list[1]}

    return ID_to_odb

# the OG2genes.tab file is parse into a dictionary where a key of odb ID will return the OG

def make_odb_to_OGset():
    
    odb_to_OGset = dict()    

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
        
    return odb_to_OGset

# the OG_dict is reversed such that a key of an OG will return a list of odb IDs that are in that group

def make_OG_to_odbSet(odb_to_OGset):

    OG_to_odbSet = dict()

    for odbID in odb_to_OGset:
        OG_list = odb_to_OGset[odbID]

        for ogID in OG_list:
            if ogID in OG_to_odbSet:
                OG_to_odbSet[ogID].add(odbID)
            else :
                OG_to_odbSet[ogID] = set()
                OG_to_odbSet[ogID].add(odbID)

    return OG_to_odbSet

# a list of OGs associated with an odb will be collected, the list of odbs associated with those OGs are added to a set
# so that a list of odbs is returned for a provided odb

def find_homologous_odb(query_list, ID_to_odb, odb_to_OGset, OG_to_odbSet):

    uniprot_to_homologs = dict()
    
    for uniprot_id in query_list:

        homologs = set()

        try:
            odb = ID_to_odb[uniprot_id]['odb']
        except KeyError:
#            print('Can not find',uniprot_id,'in ortho db')
            continue

        try:
            OGlist = odb_to_OGset[odb]
        except KeyError:
#            print('Can not find an OG for',odb)
            continue

        for OG in OGlist:
            odbList = OG_to_odbSet[OG]

            for related_odb in odbList:
                homologs.add(related_odb)

        uniprot_to_homologs[uniprot_id] = homologs

    return uniprot_to_homologs

# the homologs for each entry in the input file are found and the taxa from these homologs are written to an output file

def save_homologous_taxa_to_file(uniprot_to_homologs, taxa_by_uniprot_filename, masterset_taxa_filename):

    masterset = set()

    with open(taxa_by_uniprot_filename , 'w') as outputFile:
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
            
    with open(masterset_taxa_filename,'w') as outMaster:
        masterstr = ''

        for taxa in masterset:
            masterstr += taxa + ','
        masterstr = masterstr.rstrip(',')
        outMaster.write(masterstr)

