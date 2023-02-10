import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import hypergeom
from re import search
import csv

def csv2dict(file):
    d={}
    with open(file, "r") as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            key = row[0]
            val = list(row[1:])
            d[key] = val
    return d

def dict2csv(d,file):
    array = []
    for key, val in d.items():
        temp = [key]
        temp.extend(list(val))
        #trun = [i.split('/')[-1] for i in temp]
        #trun[0] = trun[0].split('-')[-1]
        array.append(temp) #change to trun


    with open(file, "w") as f:
        writer = csv.writer(f)
        writer.writerows(array)
        
def update_dict(d,df,key,val):
    for index, row in df.iterrows():
        if (row[key] in d) == False:
            d[row[key]]={row[val]}
        else:
            prev = d.get(row[key])
            prev.add(row[val])
            d[row[key]] = prev
    return d

def map_dict_vals(d,dict2bemapped):
    dictmapped = {}
    for s, unis in dict2bemapped.items():
        inputs = set(pd.Series(list(unis)).map(d))
        dictmapped[s] = inputs
    return dictmapped

def match_syn(x,table):
    u = table[table.Synonym.apply(lambda s: x in s)]['UniProtKB ID'].values
    if len(u):
        return u[0]
    else:
        return None
    
def split_mouse_MGI(x):
    if ':' in x:
        return 'MGI:'+ x.split(':')[1]
    else:
        return 'N.A.'
    
def split_mouse_symbol(x):
    if ':' in x:
        return x.split(':')[2]
    else:
        return 'N.A.'
    
def convert_IDs(genes, input_type):
    genes.g = genes.g.str.upper()
        
    file = '../data/simplemine_results.txt' #default is human, mouse is an option, other species not supported
    table = pd.read_csv(file,sep='\t', header=3)
    
    table['MGI'] = table['Mouse Ortholog'].apply(lambda x: split_mouse_MGI(x))
    table['Mouse Symbol'] = table['Mouse Ortholog'].apply(lambda x: split_mouse_symbol(x))
    
    table['Synonym'] = table['Synonym'].apply(lambda x: x.split(' | '))
    table['UniProtKB ID'] = table['UniProtKB ID'].apply(lambda x: x.split(' | '))
    d= pd.Series(table['UniProtKB ID'].values,index=table[input_type]).to_dict()
    genes['uniprot']=genes['g'].apply(lambda x: d.get(x))
    t = genes[genes.uniprot.isna()].copy()
    t['uniprot'] = t.g.apply(lambda x: match_syn(x,table))
    not_converted = t[t.uniprot.isna()]
    genes = genes.dropna()
    genes = pd.concat([genes,t[t.uniprot.notna()]])
    temp = genes.explode('uniprot')
    uniprot2input = pd.Series(temp.g.values, index=temp.uniprot).to_dict()
    uniprot_list = list(temp.uniprot.values)
    return uniprot_list,uniprot2input,not_converted#list(itertools.chain.from_iterable(genes.uniprot.values)), not_converted

    