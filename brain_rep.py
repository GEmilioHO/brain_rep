#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 16:49:02 2020

@author: gabrielemilioherreraoropeza
"""

# Usage: brain_rep.py [-h] -s SOURCECELL -d DESIREDCELL [-o ORGANISM]
# Example: python3 brain_rep.py -s "Neuron 001" -d "Neuron 051"

import sys, os, errno, argparse, json, shutil
import pandas as pd
import scanpy as sc

### --- Argument parser and verification

currentPATH = os.getcwd().replace('\\','/')

organisms = ['homo sapiens',
             'mus musculus'
             ]

def isFile(string):
    if os.path.isfile(string) == False:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), string)
    else:
        pass

def isCell(string):
    with open('../data/brain_cells_data.json') as json_file: 
        data = json.load(json_file)
    temp = None
    if not string in data:
        for key in data:
            if data[key]['alternative_name'] == string:
                temp = data[key]['alternative_name']
                string = key
        if temp == None:
            msg = "The cell type %s does not exist" % string
            raise argparse.ArgumentTypeError(msg)
    return string
    
def isOrg(string):
    if not string in organisms:
        msg = "The organism %s is currently not supported" % string
        raise ValueError(msg)
    return string

def isFloat(string):
	try:
		float(string)
	except:
		msg = "%s is not recognized as Float" % string
		raise argparse.ArgumentTypeError(msg)
	return float(string)
            
isFile('../data/brain_cells_data.json') # Check whether 'brain_cells_data.json' is in the correct folder
isFile('../data/gene_expression_db.csv') # Check whether 'gene_expression_db.csv' is in the correct folder

parser = argparse.ArgumentParser(
description = 'BrainRep v0.1 - a script to produce a list of predicted transcription factors to induce cellular reprogramming among brain cells',
epilog = '-s and -d are required, -o is optional as default is "homo sapiens"')

parser.add_argument('-s', 
                    '--sourceCell', 
                    dest = 'sourceCell', 
                    help = 'Source cell type name according to gene expression, see "cell_types.csv"', 
                    required = True, 
                    type = isCell
                    )

parser.add_argument('-d', 
                    '--desiredCell', 
                    dest = 'desiredCell', 
                    help='Desired cell type name according to gene expression, see "cell_types.csv"', 
                    required = True, 
                    type = isCell
                    )

parser.add_argument('-o', 
                    '--organism', 
                    dest = 'organism', 
                    default = 'homo sapiens', 
                    help = 'Organism for prediction (homo sapiens)', 
                    type = isOrg)

parser.add_argument('-t', 
                    '--levelSignificance', 
                    dest = 'levelSignificance', 
                    default = 0.05, 
                    help = 'Level of significance for t-test (0.05)', 
                    type = isFloat)

args = vars(parser.parse_args())

sourceCell = args['sourceCell']
desiredCell = args['desiredCell']
organism = args['organism']
pval_th = args['levelSignificance']

print('\nReading gene expression data...')
adata = sc.read('../data/brain_qc.h5ad') # Read gene expression data set


### --- Get scripts

import br_utils as bu


### --- Prepare out folder

if os.path.exists('../out'):
    q = input('\nDelete /out folder? (y/n/exit): ')
    if q.lower() == 'y':
        shutil.rmtree('../out')
        os.makedirs('../out')
    elif q.lower() == 'n':
        print('\nWARNING: Will not delete /out folder, that may overwrite an existing file\n')
    elif q.lower() == 'exit':
        print('\nTerminated')
        sys.exit()
        
if not os.path.exists('../out'):
    os.makedirs('../out')
    print('\nCreated folder: /out')


### --- START ---

# Convert anndata file into DataFrame
expressionDB = pd.DataFrame(data = adata.X, 
                            index = adata.obs_names,
                            columns = adata.var_names
                            )

# Print preview of expression DataFrame
print(expressionDB.head(10))

# Print conversion 
print(f'\nWill convert from {sourceCell} to {desiredCell}')

# Search for transcription-factor-coding genes in the organism
tf_list = bu.search_tf(organism = organism)

# Convert to GENENAME from the GI numbers obtained in tf_list
tf_list = bu.fetch_tf(tf_list)

# Create dataframe with expression data of only source cells
sourceCells_df = expressionDB[adata.obs.cell_type_designation_label == sourceCell]
print(f'\nExpression data from {len(sourceCells_df)} {sourceCell}s will be used...')

# Create dataframe with expression data of only desired cells
desiredCells_df = expressionDB[adata.obs.cell_type_designation_label == desiredCell]
print(f'...and expression data from {len(desiredCells_df)} {desiredCell}s will be used.')

# Find background cells for source cell type and create DataFrame
lst_bkgd_cells4sourceCells = bu.find_bkgd_cells(cellInterest = sourceCell)
bkgd4sourceCells_df = expressionDB[adata.obs.cell_type_designation_label.isin(lst_bkgd_cells4sourceCells)]

# Find background cells for desired cell type and create DataFrame
lst_bkgd_cells4desiredCells = bu.find_bkgd_cells(cellInterest = desiredCell)
bkgd4desiredCells_df = expressionDB[adata.obs.cell_type_designation_label.isin(lst_bkgd_cells4desiredCells)]

# Find differentially expressed genes
dct_diff_exp_genes_sourceCell, dct_diff_exp_genes_desiredCell = bu.find_diff_exp_genes(
    genes = adata.var_names,
    df1 = sourceCells_df,
    df2 = desiredCells_df,
    df_bkgd1 = bkgd4sourceCells_df, 
    df_bkgd2 = bkgd4desiredCells_df)

# Find which of the differentially expressed genes encode for transcription factors
dct_diff_exp_tf_sourceCell = bu.isTF(dct_genes = dct_diff_exp_genes_sourceCell, 
                                     lst_tf = tf_list)

dct_diff_exp_tf_desiredCell = bu.isTF(dct_genes = dct_diff_exp_genes_desiredCell,
                                      lst_tf = tf_list)

# Calculate gene score for all the genes of the source cell type
dct_gene_score_sourceCell = bu.get_gene_score(dct_genes = dct_diff_exp_genes_sourceCell, 
                                              cell_type = sourceCell)

# Calculate gene score for all the genes of the desired cell type
dct_gene_score_desiredCell = bu.get_gene_score(dct_genes = dct_diff_exp_genes_desiredCell,
                                               cell_type = desiredCell)

# Create a list of all the transcription factors without duplicates
combined_list_tf = list(dct_diff_exp_tf_sourceCell.keys())
combined_list_tf.extend(x for x in dct_diff_exp_tf_desiredCell.keys() if x not in combined_list_tf)

# Get interaction networks of the transcription factors
dct_interaction_network = bu.get_interaction_network(genes = combined_list_tf)

# Calculate influence scores of the transcription factors of the source cell
print('\nCalculating influence score of the transcription factors on the source cell type...')
tf_influence_score_sourceCell = bu.tf_influence_score(dct_gene_score = dct_gene_score_sourceCell, 
                                                      dct_gene_interactions = dct_interaction_network, 
                                                      dct_tf = dct_diff_exp_tf_sourceCell)

# Calculate influence scores of the transcription factors of the desired cell
print('\nCalculating influence score of the transcription factors on the desired cell type...')
tf_influence_score_desiredCell = bu.tf_influence_score(dct_gene_score = dct_gene_score_desiredCell,
                                                       dct_gene_interactions = dct_interaction_network,
                                                       dct_tf = dct_diff_exp_tf_desiredCell)

# Generate ranking of transcription factors for the source cell
print('\nGenerating ranking of transcription factors...')
tf_ranking_sourceCell = bu.rank_tf(dct_gene_score = dct_gene_score_sourceCell, 
                                   dct_tf_influence_score = tf_influence_score_sourceCell)

# Generate ranking of transcription factors for the desired cell
tf_ranking_desiredCell = bu.rank_tf(dct_gene_score = dct_gene_score_desiredCell,
                                    dct_tf_influence_score = tf_influence_score_desiredCell)

# Predicting transcription factors for conversion from source cell to desired cell
print('\nPredicting transcription factors for conversion...')
lst_predicted_tf = bu.predict_tf(dct_tf_rank_sourceCell = tf_ranking_sourceCell, 
                                 dct_tf_rank_desiredCell = tf_ranking_desiredCell, 
                                 sourceCell = sourceCell, 
                                 desiredCell = desiredCell, 
                                 df_sourceCell = sourceCells_df, 
                                 df_desiredCell = desiredCells_df,
                                 pval = pval_th)
print(f'...a total of {len(lst_predicted_tf)} transcription factors predicted.')

# Remove redundant transcription factors
print('\nRemoving redundant transcription factors...')
lst_nonredundant_tf = bu.rm_redundant_tf(dct_interaction_network = dct_interaction_network, 
                                         lst_predicted_tf = lst_predicted_tf)
print(f'...{len(lst_nonredundant_tf)} non redundant transcription factors found.')

# Read log-transformed expression data file
log_expression_df = pd.read_csv('../data/trimmed_means.csv', 
                                index_col = 'feature')

# Change cell name to alternative name
with open('../data/brain_cells_data.json') as json_file: 
    cells_file = json.load(json_file)

alt_sourceCell = cells_file[sourceCell]['alternative_name']
alt_desiredCell = cells_file[desiredCell]['alternative_name']

# Save data in .TXT file
print('\nDONE... Writing out file!')
with open(f'../out/{sourceCell}_{desiredCell}_{organism}.txt', 'w') as fileout:
    fileout.write('\t'.join(('transcription_factor', 'non-redundant_rank', 'gene_score', 'influence score', 'cpm_sourceCell', 'cpm_desiredCell')))
    fileout.write('\n')
    for n, tf in enumerate(lst_nonredundant_tf):
        fileout.write('\t'.join((tf, 
                                 str(n+1), 
                                 str(round(dct_gene_score_desiredCell[tf], 3)), 
                                 str(round(tf_influence_score_desiredCell[tf], 3)), 
                                 str(round(2**(log_expression_df[alt_sourceCell][tf]), 3)), 
                                 str(round(2**(log_expression_df[alt_desiredCell][tf]), 3)))))
        fileout.write('\n')


