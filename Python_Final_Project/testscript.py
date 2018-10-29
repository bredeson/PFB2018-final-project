#!/usr/bin/env python3

import sys
from subprocess import call
from uniprot_to_taxa_functions import *

query_file_name = sys.argv[1]
taxa_by_uniprot_filename = sys.argv[2]
masterset_taxa_filename = sys.argv[3]

query_list = query_list_parser(query_file_name)
call('date')
print('query file parsed')

ID_to_odb = make_ID_to_odb()
call('date')
print('ID_to_odb made')

odb_to_OGset = make_odb_to_OGset()
call('date')
print('odb_to_OGset made')

OG_to_odbSet = make_OG_to_odbSet(odb_to_OGset)
call('date')
print('OG_to_odbset made')

uniprot_to_homologs = find_homologous_odb(query_list, ID_to_odb, odb_to_OGset, OG_to_odbSet)
call('date')
print('uniprot to homologs are made')

save_homologous_taxa_to_file(uniprot_to_homologs, taxa_by_uniprot_filename, masterset_taxa_filename)



