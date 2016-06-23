import sys
import csv
import ipdb
#from collections import OrderedDict
import pandas as pd

i = 0
filename = []
filename_full = []

for file in sys.argv[1::]:

    filename.append(file.split('.',1)[0])
    filename_full.append(file)
    print filename_full[i]
    i += 1
    
file1 = pd.read_csv(filename_full[0], sep = '\t', usecols = [1, 4, 7], header = 0)
file2 = pd.read_csv(filename_full[1], sep = '\t', usecols = [1, 4, 7], header = 0)
file3 = pd.read_csv(filename_full[2], sep = '\t', usecols = [1, 4, 7], header = 0)

print file2.dtypes

#ipdb.set_trace()
file1['Source'] = "Consensus"
file2['Source'] = "EC"
file3['Source'] = "Raw"

merged = pd.concat([file1,file2,file3], axis=0, join='outer', join_axes=None, ignore_index=True)




#merged = file1.merge(file2, on = 'POS', how = 'outer').fillna('.')
#merged.to_csv('merged.csv', index = False)
