#!/usr/bin/env python3

import operator

overall_fq_dict = dict()

with open('GObackgroundfreq2.txt', 'r') as fileObj:
    for line in fileObj:
        term, fq = line.split('\t')
        overall_fq_dict[term] = int(fq)

#print(overall_fq_dict)

sorted_list = sorted(overall_fq_dict.items(), key=operator.itemgetter(1))

top_50 = sorted_list[-50::]
top_50_terms = []
for i in range(len(top_50)):
    top_50_terms.append[top_50[i][0]
print(top_50_terms)
                        
top_20 = sorted_list[-20::]

overall_count_dict = dict()

for key in overall_fq_dict:
    count = overall_fq_dict[key]

    if count in overall_count_dict:
        overall_count_dict[count] += 1
    else:
        overall_count_dict[count] = 1

counter = 0
for key in sorted(overall_count_dict):
    counter += 1
    #print(str(key),str(overall_count_dict[key]))

#print(sorted_list[-20::])
#print(sorted_list[-1][0])
