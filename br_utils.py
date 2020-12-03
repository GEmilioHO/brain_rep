#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 13:12:26 2020

@author: gabrielemilioherreraoropeza
"""

from Bio import Entrez
from tqdm import tqdm
from scipy.stats import ttest_ind
import os, json, requests, operator
import pandas as pd
from math import log10


### --- Identify differentially expressed genes

def find_diff_exp_genes(genes, df1, df2, df_bkgd1, df_bkgd2):
    
    print('\nFinding genes differentially expressed between the cells and the background...')
    
    ### --- Create a list for appending genes which expression is statistically significant
    dct_diff_exp_genes_df1 = {}
    dct_diff_exp_genes_df2 = {}
    
    ### --- Parse genes and identify those diff expressed between the DataFrames
    for gene in tqdm(genes):
        
        # The cell type DataFrame column corresponding to the gene that will be analysed is
        # selected
        df1_gene = df1[gene]
        df2_gene = df2[gene]
        df1_bkgd_gene = df_bkgd1[gene]
        df2_bkgd_gene = df_bkgd2[gene]
        
        # Student t-test is used to idenfity if the cells diff express a particular gene
        __, p_val_1 = ttest_ind(df1_gene,
                                df1_bkgd_gene,
                                equal_var = False,
                                nan_policy = 'omit')

        __, p_val_2 = ttest_ind(df2_gene,
                                df2_bkgd_gene,
                                equal_var = False,
                                nan_policy = 'omit')
        
        # Python rounds figures < 1e-323 to 0.0, thus in order to perform the log function
        # later on the figure all quantitites < 1e-323 are expressed as 1e-323.
        if p_val_1 == 0.0:
            p_val_1 = 1e-323
            
        if p_val_2 == 0.0:
            p_val_2 = 1e-323
         
        dct_diff_exp_genes_df1[gene] = p_val_1  
        dct_diff_exp_genes_df2[gene] = p_val_2
            
                
    #print(f'\nA total of {len(dct_diff_exp_genes_df1)} differentially expressed genes found between the source cell and the background...')
    #print(f'...while {len(dct_diff_exp_genes_df2)} differentially expressed genes between the desired cell and the background.')
    
    return dct_diff_exp_genes_df1, dct_diff_exp_genes_df2


### --- Identify transcription factors

def isTF(dct_genes, lst_tf):
    
    ### --- Create a list for appending differentially expressed transcription factors
    dct_diff_exp_tf = {}
    
    ### --- Parse genes to check whether they are transcription factors
    for key in dct_genes:
        if key in lst_tf:
            dct_diff_exp_tf[key] = dct_genes[key] 
    
    #print(f'of which {len(dct_diff_exp_tf)} are transcription factors')
    
    return dct_diff_exp_tf


### --- Search for transcription factors in an organism

def search_tf(organism):
    ### --- Insert email for accessing NCBI API
    Entrez.email = "gabriel.herrera_oropeza@kcl.ac.uk"
    
    ### --- Search all transcription factor - coding genes and obtain GI number
    # 'homo sapiens' is the default organism
    
    print('\nSearching for transcription factor - coding genes...')
    
    handle = Entrez.esearch(db = "gene", 
                            retmax = 10000, 
                            term = f"{organism}[ORGN] transcription factor"
                            )
    
    record = Entrez.read(handle)
    handle.close()
    
    return record['IdList']


### --- Convert to GENENAME from the GI numbers obtained

def fetch_tf(lst_tf):
    
    print('\nFetching genes...')
    new_lst = []
    
    ### --- Check whether a file with the conversion info exists
    # If NOT, then each GI number is fetched and converted individually
    if os.path.isfile('../data/GInum2SYMBOL.csv') == True:
        temp_dct = {}
        
        with open('../data/GInum2SYMBOL.csv', 'r') as fileopen:
            for line in fileopen:
                line = line.strip('\n')
                line = line.split(',')
                temp_dct[line[0]] = line[1]
        
        for ID in tqdm(lst_tf):
            if ID in temp_dct:
                new_lst.append(temp_dct[ID])
            else:
                handle = Entrez.efetch(db="gene", 
                                       id=ID, 
                                       rettype="gene_table", 
                                       retmode="text")
                new_lst.append(str(handle.readline().strip('\n').split(' ')[0]))
                handle.close()   
    
    if os.path.isfile('../data/GInum2SYMBOL.csv') == False:
        ### --- TODO: Parse by batched not individually
        for ID in tqdm(lst_tf):
            handle = Entrez.efetch(db="gene", 
                                   id=ID, 
                                   rettype="gene_table", 
                                   retmode="text")
            new_lst.append(str(handle.readline().strip('\n').split(' ')[0]))
            handle.close()
    
    ### --- Print total number of transcription factor - coding genes found
    print(f"\nA total of {len(new_lst)} transcription factor - coding genes found")
    
    return new_lst


### --- Find background cells for cell of interest

def find_bkgd_cells(cellInterest):
    
    with open('../data/brain_cells_data.json') as json_file: 
        cells_file = json.load(json_file)
        
    ### --- Create a list for appending identified background cells
    lst_bkgd_cells = []
    
    ### --- From alternative name get marker genes
    alt_name = cells_file[cellInterest]['alternative_name']
    alt_name = alt_name.split(' ')
    marker_gene_1 = alt_name[2]
    marker_gene_2 = alt_name[3]
    
    ### --- Get cell type
    cell_type = cells_file[cellInterest]['cell_type']
    
    ### --- Parse cells to find background cells
    for cell in cells_file:
        if cell != cellInterest:
            if not marker_gene_1 in cells_file[cell]['alternative_name'] and not marker_gene_2 in cells_file[cell]['alternative_name']:
                lst_bkgd_cells.append(cell)
            elif cell_type != cells_file[cell]['cell_type']:
                lst_bkgd_cells.append(cell)
    
    return lst_bkgd_cells


### --- Calculate gene score for the genes

def get_gene_score(dct_genes, cell_type):
    
    with open('../data/brain_cells_data.json') as json_file: 
        cells_file = json.load(json_file)
    
    cell_type = cells_file[cell_type]['alternative_name']
    
    ### --- Create dictionary for saving transcription factors as keys and gene scores as values
    dct_gene_score = {}
    
    ### --- Read .CSV file containing log-transformed fold changes as DataFrame
    df = pd.read_csv('../data/trimmed_means.csv', 
                     index_col='feature')
    
    for gene in dct_genes:
        gene_score = abs(float(df[cell_type][gene])) * (-log10(float(dct_genes[gene])))
        dct_gene_score[gene] = gene_score
        
    return dct_gene_score


### --- Get interaction network for all the transcription factors

def get_interaction_network(genes, species = 'homo sapiens'):
    
    print('\nGenerating interaction network of the transcription factors...\n')
    
    dct_species = {'homo sapiens': 9606}
    
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    # Construct URL  
    request_url = "/".join([string_api_url, output_format, method])
    
    dct_gene_interactions = {}
    
    for gene in tqdm(genes):
        
        dct_gene_interactions[gene] = {}
            
        # Set parameters      
        params = {
        
            "identifiers" : gene, # your protein
            "species" : dct_species[species], # species NCBI identifier 
            "caller_identity" : "brain_rep.py", # your app name
            "add_nodes" : 30
        
        }
        
        # Call STRING
        response = requests.post(request_url, data=params)
        
        temp_dct = {}
        
        if response.text != '':
            
            for line in response.text.strip().split("\n"):
            
                l = line.strip().split("\t")
                
                if l[0] != 'not found' and l[0] != 'Error' and len(l) > 5:
                    p1, p2 = l[2], l[3]
                    
                    # Filter the interaction according to combined score
                    combined_score = float(l[5])
                    
                    if combined_score > 0.4:
                        
                        if not p1 in temp_dct:
                            temp_dct[p1] = []
                        
                        if not p2 in temp_dct[p1]:
                            temp_dct[p1].append(p2)
                
                elif l[0] == 'not found' and l[0] == 'Error':
                    break
            
            if gene in temp_dct:
                for gene_int_1 in temp_dct[gene]:
                    dct_gene_interactions[gene][gene_int_1] = {}
                    if gene_int_1 in temp_dct:
                        for gene_int_2 in temp_dct[gene_int_1]:
                              dct_gene_interactions[gene][gene_int_1][gene_int_2] = []
                              if gene_int_2 in temp_dct:
                                  for gene_int_3 in temp_dct[gene_int_2]:
                                      dct_gene_interactions[gene][gene_int_1][gene_int_2].append(gene_int_3)
        
    return dct_gene_interactions


### --- Create function for calculating transcription factor influence score

def tf_inf_sc_func(gene_score, level, out_degree):
     
    tf_score = gene_score * (1/level) * (1/out_degree)
    
    return tf_score


### --- Calculate trancription factor influence score

def tf_influence_score(dct_gene_score, dct_gene_interactions, dct_tf):
    
    # Create dictionary to save the influence scores of the transcription factors
    dct_tf_influence_score = {}
    
    # Parse dictionary of transcription factors to calculate the influence score of each transcription factor
    for tf in tqdm(dct_tf):
        
        dct_tf_influence_score[tf] = 0
        
        for gene_1 in dct_gene_interactions[tf]:
            
            try:
                dct_tf_influence_score[tf] += tf_inf_sc_func(gene_score = dct_gene_score[gene_1],
                                                             level = 1,
                                                             out_degree = len(dct_gene_interactions[tf]))
            except KeyError:
                pass
            
            for gene_2 in dct_gene_interactions[tf][gene_1]:
                
                try:
                    dct_tf_influence_score[tf] += tf_inf_sc_func(gene_score = dct_gene_score[gene_2],
                                                                 level = 2,
                                                                 out_degree = len(dct_gene_interactions[tf][gene_1]))
                except KeyError:
                    pass
                
                for gene_3 in dct_gene_interactions[tf][gene_1][gene_2]:
                    
                    try:
                        dct_tf_influence_score[tf] += tf_inf_sc_func(gene_score = dct_gene_score[gene_3],
                                                                     level = 3,
                                                                     out_degree = len(dct_gene_interactions[tf][gene_1][gene_2]))
                    except KeyError:
                        pass
    
    return dct_tf_influence_score


### --- Get the ranking of the transcription factors

def rank_tf(dct_gene_score, dct_tf_influence_score):
    
    # Create a dictionary for saving the total score of each transcription factor
    dct_tf_rank = {}
    
    # Parse the transcription factors and calculate their total score
    for tf in dct_tf_influence_score:
        
        dct_tf_rank[tf] = dct_gene_score[tf] + dct_tf_influence_score[tf]
    
    # Sort the dictionary according to the score obtained
    SORT = sorted(dct_tf_rank.items(), key = operator.itemgetter(1), reverse = True)
    
    # Create a new dictionary for saving the sorted list of transcription factors "SORT"
    dct_tf_rank_sorted = {}
    
    # Add to the dictionary only the TOP100 transcription factors
    for n, group in enumerate(SORT):
        
        if n < 100:
            
            dct_tf_rank_sorted[group[0]] = group[1]
    
    return dct_tf_rank_sorted


### --- Check whether a specific gene between two cells is differentially expressed

def isDiff(df_sourceCell, df_desiredCell, pval):
    
    diff = False
    
    __, p_val = ttest_ind(df_sourceCell, 
                          df_desiredCell, 
                          equal_var = False, 
                          nan_policy = 'omit')
        
    if p_val < pval:
        diff = True
    
    return diff
    

### --- Predict transcription factors required for conversion

def predict_tf(dct_tf_rank_sourceCell, dct_tf_rank_desiredCell, sourceCell, desiredCell,
               df_sourceCell, df_desiredCell, pval):
    
    # Necessary transcription factors of the desired cell that are already expressed in the 
    # source cell
    lst_duplicates = [tf for tf in dct_tf_rank_desiredCell if tf in dct_tf_rank_sourceCell]
    
    # Open log-transformed expression file 
    df = pd.read_csv('../data/trimmed_means.csv', 
                     index_col = 'feature')
    
    # Change cell name to alternative name
    with open('../data/brain_cells_data.json') as json_file: 
        cells_file = json.load(json_file)
    
    sourceCell = cells_file[sourceCell]['alternative_name']
    desiredCell = cells_file[desiredCell]['alternative_name']
    
    # Create a list for appending not required transcription factors
    lst_not_required_tf = []
    
    # Search for not required transcription factors by checking their CPM
    # If CPM source cell > CPM desired cell, then it is not required
    # If CPM source cell <= PM desired cell, then check if their differentially expressed
    # from raw data
    for tf in lst_duplicates:
        log_sC = df[sourceCell][tf]
        cpm_sC = 2 ** (log_sC) # Transform to CPM from log
        log_dC = df[desiredCell][tf]
        cpm_dC = 2 ** (log_dC)
        
        if cpm_sC <= cpm_dC:
            diff = isDiff(df_sourceCell[tf], df_desiredCell[tf], pval)
            if diff == False:
                lst_not_required_tf.append(tf)
        
        elif cpm_sC > cpm_dC:
            lst_not_required_tf.append(tf)
        
    # Create a list with the predicted transcription factors
    lst_predicted_tf = [tf for tf in dct_tf_rank_desiredCell if not tf in lst_not_required_tf]
    
    return lst_predicted_tf


### --- Function for extracting all elements in a nested dictionary

def extract(dict_in):
     
    lst_out = []
    
    for key, value in dict_in.items():
        
        if isinstance(value, dict): # If value itself is a dictionary
            if not value in lst_out:
                lst_out.append(key)
            extract(value)
            
        elif isinstance(value, list): # If value is a list
            for l in value:
                if not l in lst_out:
                    lst_out.append(l)
                    
    return lst_out


### --- Determine redundant transcription factors and remove them

def rm_redundant_tf(dct_interaction_network, lst_predicted_tf):
    
    # Create dictionary for saving the transcription factors and their downstream genes
    dct_downstream_genes = {}
    
    for tf in lst_predicted_tf:
        dct_downstream_genes[tf] = extract(dict_in = dct_interaction_network[tf])
    
    # Create list for appending non-redundant transcription factors
    lst_nonredundant_tf = []
    
    for tf in lst_predicted_tf:
        temp_lst = lst_nonredundant_tf.copy()
        if len(temp_lst) > 0:
            for l in temp_lst:
                if l != tf and not tf in lst_nonredundant_tf:
                    test_list1 = dct_downstream_genes[tf]
                    test_list2 = dct_downstream_genes[l]
                    if len(test_list1) > 0:
                        
                        inc = (len([x for x in test_list1 if x in test_list2]) / len(test_list1)) * 100
                        
                        if inc < 95: # if inclusion percentage with a higher rank gene > 95 then do not append
                            lst_nonredundant_tf.append(tf)
        elif len(temp_lst) == 0:
            if len(dct_downstream_genes[tf]) > 0:
                lst_nonredundant_tf.append(tf)
    
    return lst_nonredundant_tf
    
