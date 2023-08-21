import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import hypergeom
import sys
sys.path.append('../GOCAM_Project/dev')
import os
import tqdm

import utils
import ncHGT as noncentralHGT
pd.options.display.max_colwidth = 100

def get_sizes (data): #data= dataframe with gocam IDs and gene identifiers as columns
    return data['gocam'].value_counts()
    
def get_sets (gene_list):
    sets = []
    not_in_a_set = []
    members2setID = utils.csv2dict('../data/members2setID.csv')
    setID2members_input = {}
    for g in gene_list:
        s = members2setID.get(g)
        if s != None:
            sets = sets +s
            for i in s:
                if (i in setID2members_input) == False:
                    setID2members_input[i]={g}
                else:
                    prev = setID2members_input.get(i)
                    prev.add(g)
                    setID2members_input[i] = prev
        else:
            not_in_a_set.append(g)
    return not_in_a_set, list(set(sets)),setID2members_input #remove duplicates

def filter_gene_list(gene_list, Dict):
    filtered_gene_list = []
    filtered_out = []
    for gene in gene_list:
        if gene in Dict:
            filtered_gene_list.append(gene)
        else:
            filtered_out.append(gene)
    return filtered_out, filtered_gene_list

def count_genes(gene_list, Dict):
    gocam_counts = {} #key=gocam, value=list of genes in gocam that are also in the user's list
    for g in gene_list: 
            gocams = Dict.get(g)
            for gocam in gocams:
                if (gocam in gocam_counts) == False:
                    gocam_counts[gocam]=[g]
                else:
                    prev = gocam_counts.get(gocam)
                    prev.append(g)
                    gocam_counts[gocam] = prev
    return gocam_counts

#BENJAMINI HOCHBERG CORRECTION applied in correct_pval_and_format()
#ncHGT is either False (indicating that regular HGT should be done) or a positive integer denoting N for ncHGT
def hgt(counts, gocam_sizes, alpha, gene_list_size, background_gene_list_size, ncHGT = False):
    """ ncHGT is either False (for set or standard methods) or corresponds to N """
    results = []
    iterator = tqdm.tqdm(counts.items())
    for gocam, gene_list in iterator:
        count = len(gene_list) 
        gocam_size = gocam_sizes[gocam]
        pvalue = None
        if ncHGT:
            if count <=1: #avoid unnecessary calls to BiasedUrn due to computation time
                pvalue = 1
            else:
                pvalue = noncentralHGT.do_ncHGT(count -1,gocam,background_gene_list_size,ncHGT)
        else: #set or standard methods
            pvalue = hypergeom.sf(count-1, background_gene_list_size,  gocam_size, gene_list_size) 
        if pvalue < 1: #alpha:
            r = (gocam, pvalue, count, gocam_size, gene_list )
            results.append(r)
    return results

#Benjamini Hochberg correction
def correct_pval_and_format(enriched_gocams, background_num_gocams,show_significant,alpha):
    df = pd.DataFrame(enriched_gocams, columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df.sort_values('pval (uncorrected)',inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['FDR_val'] = (df.index+1)*alpha/background_num_gocams
    df['Less_than'] = (df['pval (uncorrected)'] < df['FDR_val'])
    index = df.Less_than.where(df.Less_than==True).last_valid_index()
    df_significant = df
    if (show_significant):
        df_significant = df.loc[0:index].copy()
        if index == None:
            df_significant = pd.DataFrame(columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df_display = df_significant[['url','pval (uncorrected)', '# entities in list', '#entities in model','shared entities in gocam']].copy()
    #modelID2title = pd.read_csv('../data/modelID2title_mouse.csv')
    temp = pd.read_csv('../data/modelID2title_mouse.csv',header = 0,names=['gocam','title'])
    modelID2title = pd.Series(temp.title.values,index=temp.gocam).to_dict()
    df_display['title'] = df_display['url'].map(modelID2title)
    cols = df_display.columns.to_list()
    cols[0]='title'
    cols[-1]='url'
    df_display = df_display[cols]
    return df_display

#Dict can only contain 1 instance of each gene per gocam (no duplicates)
#show_significant only affects the multiple testing correction. If the uncorrected pval > alpha, hgt() will already remove it
def enrich(gene_list, uni_list,uniprot2input,gocam_sizes, Dict, ncHGT=False, show_significant=True,alpha=.05):
    background_gene_list_size = len(Dict)
    if ncHGT: 
    #we consider the background size to be equal to the total # of genes 
    #(the sum of the weights of all entities would double count genes that occur in multiple sets
    #... is this the right thing to do though?
        background_gene_list_size = len(utils.csv2dict('../data/ID2gocam_mouse_ff.csv'))
        
    not_in_a_set, sets, setID2members_input_uni = get_sets(uni_list)
    
    setID2members_input = utils.map_dict_vals(uniprot2input, setID2members_input_uni)
    
    filtered_out1, set_list_filtered = filter_gene_list(sets,Dict)
    filtered_out2, gene_list_filtered = filter_gene_list(uni_list, Dict) #need to clean gene_list to only include genes in the gocam
    
    
    filtered_list = gene_list_filtered + set_list_filtered
    gene_list_size = len(filtered_list)
    
    flist2input = {**uniprot2input, **setID2members_input}
    filtered_list_as_genes = set(pd.Series(list(filtered_list)).map(flist2input).explode())
    filtered_out_genes = set(gene_list) - filtered_list_as_genes
    
    counts = count_genes(filtered_list, Dict)
    
    N_ncHGT = False
    if ncHGT == True:
        N_ncHGT = len(gene_list)-len(filtered_out_genes)
        if N_ncHGT <= 0:
            return "error no genes found in gocams"
        
    enriched_gocams = hgt(counts, gocam_sizes, alpha, gene_list_size, background_gene_list_size, ncHGT=N_ncHGT)
    background_num_gocams = len(gocam_sizes)
    df_display = correct_pval_and_format(enriched_gocams, background_num_gocams,show_significant,alpha)
    return filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display
    
def enrich_wrapper(filename, id_type, return_all = False, method = 'set', show_significant=True,alpha=.05,fpath= '../test_data', display_gene_symbol = True):
    """ returns (gene_list, filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)"""
    #set method files
    gcs = '../data/gocam_sizes_mouse.csv'
    id2g = '../data/ID2gocam_mouse.csv'
    
    #standard method files
    if method == 'standard':
        gcs = '../data/gocam_sizes_mouse_ff.csv'
        id2g = '../data/ID2gocam_mouse_ff.csv'
    
    gene_list = pd.read_csv(os.path.join(fpath,filename),header=None,names = ['g'])
    
    #normally not needed, but I found a bug where HSPA1A and HSPA1B are listed as synonyms, both in Simplemine and official sources like the Alliance
    gene_list.drop_duplicates(inplace = True) 
    
    gene_list_converted = []
    uniprot2input = {}
    not_converted = []
    
    #conversion to uniprot IDs not needed for a list of uniprot IDs
    if id_type == 'uniprot':
        gene_list_converted = gene_list.g
        uniprot2input = pd.Series(gene_list_converted.values,index=gene_list_converted).to_dict()
    else:
        gene_list_converted, uniprot2input, not_converted = utils.convert_IDs(gene_list,id_type)
    
    #read in dictionary and the gocam sizes
    x = pd.read_csv(gcs)
    gocam_sizes = pd.Series(x.sizes.values,index=x.gocam)
    Dict = utils.csv2dict(id2g)
    
    #call enrich()
    ncHGT = False
    if method == 'ncHGT':
        ncHGT = True
    #results: (filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
    results = enrich(list(gene_list.g), gene_list_converted, uniprot2input, gocam_sizes, Dict, ncHGT = ncHGT, show_significant = show_significant, alpha=alpha)
    
    if display_gene_symbol == True:
        results[4]['shared entities in gocam'] = utils.uniprot2gene(results[4]['shared entities in gocam'])
    if method == 'set' or method == 'ncHGT':
        print(f"Analysis run on {len(results[1])} entities from {len(gene_list)-len(results[0])} out of {len(gene_list)} input genes")
    if return_all:
        return (gene_list, *results)
    else:
        return results[4]

        
def compare2standard(d1, filename , symbol, alpha=.05,fpath = '../test_data/processed'):
    
    d2= enrich_wrapper(filename,symbol,method='standard',alpha = alpha,fpath = fpath)


    r1 = set(d1.title.values)
    r2 = set(d2.title.values)

    dif12 = set(r1)-set(r2)
    dif21 = set(r2)-set(r1)
        
    d2_all = enrich_wrapper(filename,symbol,method='standard',show_significant= False, alpha = alpha,fpath = fpath)
    dif12 = list(dif12)
    sizes = {}
    for i in dif12:
        s = d2_all[d2_all.title == i]["#entities in model"].values[0]
        sizes[i] = s

    df = d1[d1.title.apply(lambda x: x in dif12)].copy()
    if len(df) > 0:
        
        df['gene list size'] = df.title.apply(lambda x: sizes.get(x))
    
    print(f"Standard method yields {len(r2)} results, {len(dif21)} of which are unique")
    
    return df
        
        