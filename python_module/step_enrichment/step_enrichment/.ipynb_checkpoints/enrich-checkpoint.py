import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import hypergeom
import os
import os.path
import tqdm

from .utils import csv2dict, map_dict_vals, convert_IDs, reverse_dict, backend2gene, get_data_directory
from .ncHGT import do_ncHGT 
pd.options.display.max_colwidth = 100


    
def get_sets(gene_list):
    """map list of genes to all sets that contain members of that list"""
    dpath = get_data_directory()

    sets = []
    not_in_a_set = []
    members2setID = csv2dict(os.path.join(dpath, 'members2setID.csv'))
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
    """remove members of gene_list that are not in Dict.
    use function to filter a user's input list of genes based on those that appear at least 
    once in the gocam model database"""
    filtered_gene_list = []
    filtered_out = []
    for gene in gene_list:
        if gene in Dict:
            filtered_gene_list.append(gene)
        else:
            filtered_out.append(gene)
    return filtered_out, filtered_gene_list

def count_genes(gene_list, Dict):
    """ count number of genes in user's gene_list that are in each gocam"""
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
def hgt(counts, gocam_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT = False, 
        kegg = False, input_type = None, backend_type = None, enrich_against = None):
    """ performs either the hypergeometric test or our introduced test using Fisher's noncentral hypergeometric dist.
    Whether our unweighted set enrichment or the standard HGT is performed is determined upstream based on what
    Dict of gocams->entities and filtered gene_list are passed into count_genes().
    ncHGT is either False (for set or standard methods) or corresponds to N """
    
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
                #changed from count -1 to count on 7/15/24. This error would not invalidate any results in the paper. It would make them more significant, as 
                #p-values are the probability of obtaining something "as extreme as" k. Survival functions give P(K > k) which is why we use hypergeom.sf(count-1) 
                #in the next code block. However, do_ncHGT computes P(K >= k). This error led to reporting values as slightly less significant by also adding 
                #P(K = k - 1) to the sum
                pvalue = do_ncHGT(count,gocam,background_gene_list_size,ncHGT, 
                                                kegg = kegg, input_type = input_type, backend_type = backend_type, enrich_against = enrich_against)
        else: #set or standard methods
            pvalue = hypergeom.sf(count-1, background_gene_list_size,  gocam_size, gene_list_size) 
        if pvalue < 1: #FDR:
            r = (gocam, pvalue, count, gocam_size, gene_list )
            results.append(r)
    return results

#Benjamini Hochberg correction
def correct_pval_and_format(enriched_gocams, background_num_gocams,FDR,kegg, enrich_against = None):
    """performs Benjamini Hochberg correction to control the false discovery rate and formats output for display"""
    dpath = get_data_directory()

    df = pd.DataFrame(enriched_gocams, columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df.sort_values('pval (uncorrected)',inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['FDR_val'] = (df.index+1)*FDR/background_num_gocams
    df['Less_than'] = (df['pval (uncorrected)'] < df['FDR_val'])
    index = df.Less_than.where(df.Less_than==True).last_valid_index()
    df_significant = df
    
    df_significant = df.loc[0:index].copy()
    if index == None:
        df_significant = pd.DataFrame(columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df_display = df_significant[['url','pval (uncorrected)', '# entities in list', '#entities in model','shared entities in gocam']].copy()
    #modelID2title = pd.read_csv('../data/modelID2title_mouse.csv')
    
    temp = pd.read_csv(os.path.join(dpath, 'modelID2title_mouse.csv'),header = 0,names=['pathway','title'])
    if kegg:
        temp = pd.read_csv(os.path.join(dpath,f'kegg/{enrich_against}_name.txt'),header = None, sep = '\t', names=['pathway','title'])
    modelID2title = pd.Series(temp.title.values,index=temp.pathway).to_dict()
    
    df_display['title'] = df_display['url'].map(modelID2title)
    cols = df_display.columns.to_list()
    cols[0]='title'
    cols[-1]='url'
    df_display = df_display[cols]
    return df_display

#Dict can only contain 1 instance of each gene per gocam (no duplicates)
def _enrich(gene_list, uni_list,uniprot2input,pathway_sizes, Dict, ncHGT=False,FDR=.05, 
           kegg = False, input_type = None, backend_type = None, enrich_against = None):
    """uni_list is the list of uniprot IDs, because the backend dictionary, Dict, is gocam_id-> list(uniprot id's).
    uniprot2input is a dictionary keeping track of which of the user's inputs mapped to which uniprot id's so results can be 
    displayed in the user's inputted format, as the mapping is not always 1:1."""
    
    dpath = get_data_directory()
    background_gene_list_size = len(Dict)
    if kegg == False:
        if ncHGT: 
        #we consider the background size to be equal to the total # of genes 
        #(the sum of the weights of all entities would double count genes that occur in multiple sets
        #... is this the right thing to do though?
            background_gene_list_size = len(csv2dict(os.path.join(dpath,'ID2gocam_mouse_ff.csv')))

        not_in_a_set, sets, setID2members_input_uni = get_sets(uni_list)

        setID2members_input = map_dict_vals(uniprot2input, setID2members_input_uni)

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

        enriched_gocams = hgt(counts, pathway_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT=N_ncHGT)
        background_num_gocams = len(pathway_sizes)
        df_display = correct_pval_and_format(enriched_gocams, background_num_gocams,FDR, kegg)
        return filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display
    
    else:
        if ncHGT: 
            background_gene_list_size = len(csv2dict(os.path.join(dpath,f'kegg/{input_type}-to-reaction.map'), sep = '\t')) #should adjust this when filtering is applied
        gene_list_size = len(uni_list) #should be named backend_list. number of reactions
        counts = count_genes(uni_list, Dict) #number of reactions per pathway
        counts = remove_non_pathways_kegg(counts, enrich_against)
        N_ncHGT = False
        if ncHGT == True:
            N_ncHGT = len(gene_list) #number of inputted genes #CHECK
        enriched_pathways = hgt(counts, pathway_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT=N_ncHGT, 
                                kegg = kegg, input_type = input_type, backend_type = backend_type, enrich_against = enrich_against)
        background_num_pathways = len(pathway_sizes)
        df_display = correct_pval_and_format(enriched_pathways, background_num_pathways,FDR, kegg, enrich_against = enrich_against)
        
        return ['no filtering with kegg'],uni_list,{'not done with kegg':'not done with kegg'}, uniprot2input, df_display
        

def construct_dicts(custom):
    pass

def remove_non_pathways_kegg(counts, enrich_against):
    """kegg has a number of pathways that aren't meaningful pathways, such as ko01100 'Metabolic Pathways'. These are also very large and slow down the code"""
    if enrich_against == 'pathways':
        for p in ['ko01100','ko01110','ko01120', 'ko01200']:
            counts.pop(p, None)
    elif enrich_against == 'modules':
        pass
    return counts


def enrich(filename, input_type, method = 'set', return_all = False, FDR=.05,fpath= '', display_gene_symbol = True, display_input = False, 
                   kegg = False, enrich_against = 'module', custom = {}):
    """ wrapper to perform enrichment given a filename, gene ID type, enrichment method, and false discovery rate.
    other parameters:
    valid gene IDs for GOCAMs/Reactome: 'Gene Symbol', 'ENSEMBL ID', 'uniprot'
    valid gene IDs for kegg: 'ko', 'enzyme' (EC number)
    if enriching with kegg, set kegg = True.
    - enrich against can be 'module' or 'pathways'
    return_all: 
        if false, only returns the dataframe displaying results. 
        if true: returns (gene_list, filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
        return_all = True is not just for debugging. User may want to know which of their input genes were filtered out as well as how
        the IDs were mapped, as uniprot IDs can sometimes map to more than one HGNC gene symbol
    display_gene_symbol: if true, display HGNC symbols on output regardless of input ID type"""
    
    backend_type = 'reaction'
    if (backend_type != 'reaction' or enrich_against != 'module') and kegg == False: #planning to allow other backend types in future. 
        raise ValueError("backend_type, enrich_against argument(s) should only be used with kegg")
    
    dpath = get_data_directory()
    #set method files
    gcs = os.path.join(dpath,'gocam_sizes_mouse.csv')
    id2g = os.path.join(dpath,'ID2gocam_mouse.csv')
    sep = ','
    #standard method files
    if method == 'standard':
        gcs = gcs[:-4] + '_ff.csv'   #'../data/gocam_sizes_mouse_ff.csv'
        id2g = id2g[:-4] + '_ff.csv' #'../data/ID2gocam_mouse_ff.csv'
        
    if kegg:
        #gcs = '../data/kegg/pathway_sizes_kegg.csv'
        print('kegg functionality: map kegg orthologs or ECs to reactions and enrich reactions against pathways or modules. Reactions are treated as "sets." Custom addition or removal of entities not implemented yet. Weighted enrichment against pathways is slow due to large sizes; enrich_against for this is set to "module" as default. Certain non-meaningful pathways such as k01100 "Metabolic Pathways" are removed from enrichment.')
        sep = '\t'
        id2g = os.path.join(dpath,f'kegg/{backend_type}-to-{enrich_against}.map')
        if method == 'standard':
            raise ValueError('standard mapping not implemented yet for kegg. use set or ncHGT')
        
        
    if custom != {}:
        id2g = construct_dicts(custom)
        
    
    
    gene_list = pd.read_csv(os.path.join(fpath,filename),header=None,names = ['g'])
    
    #normally not needed, but I found a bug where HSPA1A and HSPA1B are listed as synonyms, both in Simplemine and official sources like the Alliance
    gene_list.drop_duplicates(inplace = True) 
    
    gene_list_converted = []
    backend2input = {}
    not_converted = []
    
    #conversion to uniprot IDs not needed for a list of uniprot IDs
    #kegg does not use conversion. input must be ko or enzyme, and they are converted to reactions as an analog for sets
    if input_type in ['uniprot']:
        gene_list_converted = gene_list.g
        backend2input = pd.Series(gene_list_converted.values,index=gene_list_converted).to_dict()
    else:
        gene_list_converted, backend2input, not_converted = convert_IDs(gene_list,input_type, backend_type = backend_type, enrich_against = enrich_against)
    #read in dictionary and the gocam sizes
    #x = pd.read_csv(gcs)
    #gocam_sizes = pd.Series(x.sizes.values,index=x.gocam)
    Dict = csv2dict(id2g,sep = sep)
    pathways2members = reverse_dict(Dict)
    sizes = []
    for k,v in pathways2members.items():
        sizes.append(len(v))
    
    gocam_sizes = pd.Series(sizes, index = pathways2members.keys())    
    ncHGT = False
    if method == 'ncHGT':
        ncHGT = True
    #results: (filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
    results = _enrich(list(gene_list.g), gene_list_converted, backend2input, gocam_sizes, Dict, ncHGT = ncHGT, FDR=FDR, 
                     kegg = kegg, input_type = input_type, backend_type = backend_type,enrich_against= enrich_against)
    
    if display_gene_symbol == True:
        results[4]['shared entities in gocam'] = backend2gene(results[4]['shared entities in gocam'])
        results[4]['shared entities in gocam'] = results[4]['shared entities in gocam'].apply(lambda x: [x_.replace('sset:','set:') for x_ in x])
    if display_input == True:
        #overrides default of display_gene_symbol
        results[4]['shared entities in gocam'] = results[4]['shared entities in gocam'].apply(lambda x: [backend2input.get(x_) for x_ in x])

    if method == 'set' or method == 'ncHGT':
        print(f"Analysis run on {len(results[1])} entities from {len(gene_list)-len(results[0])} out of {len(gene_list)} input genes")
    

    if return_all:
        return (gene_list, *results)
    else:
        return results[4]
    