'''
Calculates the area-under-the-curve for user-specified gene sets in given datasets
Author: Mudra Hegde
E-mail: mhegde@broadinstitute.org
'''
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.metrics import auc
from scipy.stats import rankdata
import argparse, csv, os, warnings

warnings.filterwarnings('ignore')

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
        type=str,
        help='.txt input file with sgRNAs in first column and LFCs in following columns')
    parser.add_argument('--chip-file',
        type=str,
        help='.txt chip file to map sgRNAs to gene symbols;First column should be sgRNAs and second column should be gene symbols')
    parser.add_argument('--gene-set-folder',
        type=str,
        help='Folder with files of gene sets;One column of genes with no header')
    parser.add_argument('--outputfile',
        type=str,
        help='Outputfile')
    return parser

def calc_auc(df,column):
    df_ecdf = ECDF(df[column])
    df_x = np.linspace(min(df[column]),max(df[column]))
    df_y = df_ecdf(df_x)
    df_x = np.insert(df_x,0,0)
    df_y = np.insert(df_y,0,0)
    df_x = np.append(df_x,1.0)
    df_y = np.append(df_y,1.0)
    df_auc = auc(df_x,df_y)
    return df_auc

if __name__ == '__main__':
    args = get_parser().parse_args()
    input_df = pd.read_table(args.input_file)
    in_colnames = list(input_df.columns)
    in_colnames[0] = 'Barcode Sequence'
    input_df.columns = in_colnames
    chip_df = pd.read_table(args.chip_file)
    ch_colnames = list(chip_df.columns)
    ch_colnames[0:2] = ['Barcode Sequence','Gene Symbol']
    chip_df.columns = ch_colnames 
    merged = pd.merge(input_df,chip_df,on='Barcode Sequence')
    outputfile = args.outputfile
    gene_set_folder = args.gene_set_folder
    gene_sets = os.listdir(gene_set_folder)
    gene_set_labs = [x[:-4] for x in gene_sets]
    output_cols = ['Sample']
    output_cols.extend(gene_set_labs)
    cols = list(input_df.columns)[1:]
    with open(outputfile,'w') as o:
        w =  csv.writer(o,delimiter='\t')
        w.writerow((output_cols))
        for c in cols:
            print c
            df_c = merged[['Gene Symbol',c]]
            df_c.columns = ['Gene Symbol','sgRNA measured value']
            df_c.loc[:,'percent_rank'] = rankdata(df_c['sgRNA measured value'])/float(len(df_c))
            aucs = []
            row = [c]
            for i,g in enumerate(gene_sets):
                gene_set = pd.read_table(gene_set_folder+'/'+g,header=None)
                gene_set.columns = ['Gene Symbol']
                df_gs = pd.merge(df_c,gene_set,on='Gene Symbol')
                auc_gs = calc_auc(df_gs,'percent_rank')
                aucs.append(auc_gs)
            row.extend(aucs)
            w.writerow((row))