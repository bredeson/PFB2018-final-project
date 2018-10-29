#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt  
#import pandas as pd  
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage

total_taxa_file = sys.argv[1]
taxa_by_protein = sys.argv[2]

ID_to_taxa_dict = dict()

ID_to_taxa_array = dict()

with open(total_taxa_file, 'r') as fileObj:
    for line in fileObj:
        list_of_taxa = line.split(',')

with open(taxa_by_protein, 'r') as fileObj:
    for line in fileObj:
        id = line.split()[0]
        taxaset = line.split()[1].split(',')
        ID_to_taxa_dict[id] = taxaset

for ID in ID_to_taxa_dict:
    check_list = []

    for taxa in list_of_taxa:
        if taxa in ID_to_taxa_dict[ID]:
            check_list.append(1)
        else:
            check_list.append(0)

    ID_to_taxa_array[ID] = check_list

matrix_str = ''
for key in ID_to_taxa_array:
    array = ID_to_taxa_array[key]
    array_str = str(array)
    array_str = array_str.rstrip(']')
    array_str = array_str.lstrip('[')
    matrix_str += array_str+'; '
        
matrix_str = matrix_str.rstrip()
matrix_str = matrix_str.rstrip(';')

taxa_matrix = np.array(np.mat(matrix_str))

z = linkage(taxa_matrix,metric='cosine')

plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.savefig('my_dend.png')
   
        
