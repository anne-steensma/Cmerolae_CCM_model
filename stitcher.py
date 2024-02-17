import sys
import pandas as pd

start_filename = 'Results_1.csv'
full_dataset = pd.read_csv(start_filename,index_col=0)

for i in range(2,601):
    filename = 'Results_{0}.csv'.format(i)
    partial_dataset = pd.read_csv(filename,index_col=0)
    full_dataset = full_dataset.append(partial_dataset,ignore_index=True)

full_dataset.to_csv('Full_Dataset.csv')