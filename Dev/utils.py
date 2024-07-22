import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import hypergeom
from re import search
import csv
import os.path

def csv2dict(file, sep = ','):
    d={}
    with open(file, "r") as f:
        csv_reader = csv.reader(f, delimiter = sep)
        for row in csv_reader:
            key = row[0]            
            val = list(row[1:])
            d[key] = val
    return d

def dict2csv(d,file, sep = ','):
    array = []
    for key, val in d.items():
        temp = [key]
        temp.extend(list(val))
        #trun = [i.split('/')[-1] for i in temp]
        #trun[0] = trun[0].split('-')[-1]
        array.append(temp) #change to trun


    with open(file, "w") as f:
        writer = csv.writer(f, delimiter = sep)
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

def map_dict_vals(mapping_dict,dict2bemapped):
    dictmapped = {}
    for k, vals in dict2bemapped.items():
        mapped_vals = set(pd.Series(list(vals)).map(mapping_dict))
        dictmapped[k] = inputs
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
    
def convert_IDs(genes, input_type, backend_type = 'reaction', enrich_against = None):
    """converts input type to uniprot IDs for gocams"""
    if input_type == 'ko' or input_type == 'enzyme': #kegg
        file = f'../data/kegg/wol2/{input_type}-to-{backend_type}.map'
        if input_type == 'enzyme' and os.path.isfile(file) == False:
            kegg_make_EC_reaction()
        d = csv2dict(file, sep = '\t')
        genes[backend_type]=genes['g'].apply(lambda x: d.get(x))
        t = genes[genes[backend_type].isna()].copy()        
        not_converted = t[t[backend_type].isna()]
        genes = genes.dropna()
        genes = pd.concat([genes,t[t[backend_type].notna()]])
        temp = genes.explode(backend_type)
        temp.drop_duplicates(subset=backend_type,inplace = True)
        
        #remove reactions that are not annotated to a pathway
        d2 = csv2dict(f'../data/kegg/wol2/{backend_type}-to-{enrich_against}.map', sep = '\t')
        mask = temp[backend_type].apply(lambda x: True if d2.get(x) != None else False)
        temp = temp[mask]
        
        reactions2input = pd.Series(temp.g.values, index=temp[backend_type]).to_dict()
        reaction_list = list(temp[backend_type].values)
        return reaction_list, reactions2input, not_converted
    else:
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
        temp.drop_duplicates(subset='uniprot',inplace = True)
        uniprot2input = pd.Series(temp.g.values, index=temp.uniprot).to_dict()
        uniprot_list = list(temp.uniprot.values)
        return uniprot_list,uniprot2input,not_converted#list(itertools.chain.from_iterable(genes.uniprot.values)), not_converted

def u2ghelper(x,d):
    x_new = []
    for u in x:
        x_new.append(d.get(u,u))
    return x_new

def backend2gene(series, kegg = False):
    if kegg:
        print('not implemented yet')
    else:
        file = '../data/simplemine_results.txt' #default is human, mouse is an option, other species not supported
        table = pd.read_csv(file,sep='\t', header=3)
        table['UniProtKB ID'] = table['UniProtKB ID'].apply(lambda x: x.split(' | '))
        table = table.explode('UniProtKB ID')
        d= pd.Series(table['Gene Symbol'].values,index=table['UniProtKB ID']).to_dict()

        return series.apply(lambda x: u2ghelper(x,d))

def reverse_dict(original_dict):
    reversed_dict = {}
    
    for key, values in original_dict.items():
        for value in values:
            if value not in reversed_dict:
                reversed_dict[value] = []
            reversed_dict[value].append(key)
    
    return reversed_dict

def kegg_make_EC_reaction():
    # Read the input file
    with open('../data/kegg/wol2/reaction_enzyme.txt', 'r') as file:
        lines = file.readlines()

    # Process the lines to create the desired output format
    output_lines = []
    for line in lines:
        identifier, ec_numbers = line.strip().split('\t')
        ec_numbers = ec_numbers.strip("[]").replace("'", "").split(', ')
        output_line = identifier + '\t' + '\t'.join(ec_numbers)
        output_lines.append(output_line)

    # Write the output to a TSV file
    with open('../data/kegg/wol2/reaction-to-enzyme.map', 'w') as file:
        for output_line in output_lines:
            file.write(output_line + '\n')

    print("reaction-to-enzyme.map created from reaction_enzyme.txt")
    
    re_ec = csv2dict('../data/kegg/wol2/reaction-to-enzyme.map', sep = '\t')
    ec_re = reverse_dict(re_ec)
    dict2csv(ec_re, '../data/kegg/wol2/enzyme-to-reaction.map', sep = '\t')
    print("enzyme-to-reaction.map created")
    
def kegg_make_pathway_reaction(backend_type, enrich_against):
    
    be_pa = csv2dict(f'../data/kegg/wol2/{backend_type}-to-{enrich_against}.map', sep = '\t')
    pa_be = reverse_dict(be_pa)
    dict2csv(pa_be, f'../data/kegg/wol2/{enrich_against}-to-{backend_type}.map', sep = '\t')
    print(f"{enrich_against}-to-{backend_type}.map created")