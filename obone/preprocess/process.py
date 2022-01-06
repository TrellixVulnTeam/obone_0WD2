import GEOparse
import pandas as pd
import sqlite3
import csv
import os
import sys
import math
import re
import tempfile
import statistics
import argparse
import pymongo
import gzip
import shutil
import numpy as np
import io
import scanpy as sc
import glob
from anndata import AnnData

accessionID = ''        # GSE ID
experimentType = ''     # m or rna      (microarray or RNA-Seq)
norm = ''               # cpm/None  Need to normalize data
takeLog = False         # True/False    need to take log of data (data is not log normalized yet)
sep_rna_files = ''      # True/False    RNA-Seq data in separate files (dir of files for each sample or 1 file for all samples)?
all_files = True        # True/False    True: generate all files; False: only generate surv and ih
make_expr = True        # True/False    True: need to run rna_analyze(); False: generated log normalized expression file in Jupyter

if accessionID == '' or experimentType == '' or norm == '' or sep_rna_files == '':
    print("System usage: [GEO Accession ID] [m (microarray) or rna (RNA-Seq)] [takelog: True/False] [Normalization: cpm or tpm or none] [raw seq: True/False]")
    quit()

# GEOParse

path='./'+str(accessionID)+'_family.soft.gz'
path_dir='./'

if os.path.exists(path):
    # Load from existing file
    print("-Loading from", path)
    gse = GEOparse.get_GEO(filepath = path, silent=False)

else:
    # Download GSE from GEO and load it
    print('-Downloading', str(accessionID))

    gse = GEOparse.get_GEO(geo=str(accessionID), destdir=path_dir)
    
print("Done processing")

# Unzip SOFT file
with gzip.open(path_dir+str(accessionID)+'_family.soft.gz', 'rb') as f_in:
    with open(path_dir+str(accessionID)+'_family.soft.txt', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


"""

Generates expr.txt file for microarray experiments

"""

def m_analyze(gpl):
    print('Starting m_analyze')
    gse_df = pd.DataFrame()
    num_rows = 0
    for name, gsm in gse.gsms.items():
        if gsm.metadata['platform_id'][0] != gpl.name:
            continue
        print(name)
        column = gsm.table['VALUE']
        gse_df[name] = column
        
        if num_rows == 0:
            num_rows = len(column)
        elif num_rows != len(column):
            print('Warning: num of probes for % != other sample'%name)
        
        if takeLog == True:
            print("Taking log")
            gse_df = gse_df.apply(lambda x: np.log2(x+1))

    gse_df.insert(0, 'ID', gsm.table['ID_REF'])
    gse_df.to_csv(path_dir+str(accessionID)+'gse_withid-%s.txt'%gpl.name, header=True, index=False, sep='\t')


    if os.path.exists(path_dir+str(accessionID)+'-%s.txt'%gpl.name):
        os.remove(path_dir+str(accessionID)+'-%s.txt'%gpl.name)

    # Get gene names for GPL IDs
    print(gpl.table.columns)

    col_symbol = []
    col_title = []
    
    if 'gene_assignment'in gpl.table.columns:
        lst = gpl.table['gene_assignment']
        for l in lst:

            if l != '---':
                l = str(l)
                x = [y.split(' // ') for y in l.split(' /// ')]
                col_title.append(" /// ".join(i[2] if len(i) > 2 else "ERR" for i in x))
                col_symbol.append(" /// ".join(i[1] if len(i) > 2 else "ERR" for i in x))
            else:
                col_symbol.append("---")
                col_title.append("---")
                
    # check it, Definition and symbol can be saved with different name other than the following
    if 'GeneSymbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'GeneSymbol':'Symbol'})
    if 'Gene Symbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Symbol':'Symbol'})
    if 'gene symbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'gene symbol':'Symbol'})

    if 'geneName'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'geneName':'Definition'})
    if 'Gene Name'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Name':'Definition'})
    if 'Gene Title'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Title':'Definition'})
    if 'Gene Description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Gene Description': 'Definition'})
    if 'Description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Description': 'Definition'})
    if 'DESCRIPTION' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'DESCRIPTION': 'Definition'})
    if 'GENE_SYMBOL' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'GENE_SYMBOL': 'Symbol'})
    if 'GENE_SYM' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'GENE_SYM': 'Symbol'})
    if 'description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'description': 'Definition'})
    if 'Gene_Name' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Gene_Name': 'Symbol'})
    if 'Gene_Desc' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Gene_Desc':'Definition'})


    if 'gene_assignment'in gpl.table.columns:
        gpl.table['Symbol'] = col_symbol
        gpl.table['Definition'] = col_title
    
    gpl.table.to_csv(path_dir+str(accessionID)+'-%s.txt'%gpl.name, header=True, index=False, sep='\t')

    # Generate expr.txt
    gpl_df = pd.read_csv(str(accessionID)+'-%s.txt'%gpl.name, sep='\t')
    gpl_df = gpl_df[['ID', 'Symbol', 'Definition']]
    data_expr = gpl_df.merge(gse_df, how='right', on='ID') 
    data_expr.insert(1, 'Name', data_expr['Symbol'] + ':' + data_expr['Definition'])
    del data_expr['Symbol']
    del data_expr['Definition']

    data_expr = data_expr.rename({'ID':'ProbeID'},axis=1)

    data_expr.to_csv(path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name), header=True, index=False, sep='\t')
    print("Done with m_analyze")


"""

Generates expr file for RNA-seq experiments

"""

def rna_analyze(gpl):
    if sep_rna_files == 'True':
        files = glob.glob(path_dir+'raw/*.txt.gz')
        print('j', files)
        df = None
        for f in files:
            arr = re.sub(".*(GSM[0-9]+).*", "\\1", f)
            df1 = pd.read_csv(f, sep="\t", names=['ProbeID', arr])
            df1 = df1[~df1['ProbeID'].isin(['no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned', 'alignment_not_unique'])]
            if df is None:
                df = df1
            else:
                if arr in df.columns:
                    df[arr] += df1[arr]
                else:
                    df = df.merge(df1, 'outer', on='ProbeID')
        
        # Add gene names
        scaff_path = '/booleanfs2/sahoo/Data/SeqData/genome/'
        if df['ProbeID'].str.contains('ENSG').any():
            gene_names = pd.read_csv(scaff_path+'Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt', sep='\t', names=['T', 'ProbeID', 'Name'])
            del gene_names['T']
            df['ProbeID'] = df['ProbeID'].str.split('.', expand=True)[0]
            df = df.merge(gene_names, how='left', on='ProbeID')

        elif df['ProbeID'].str.contains('ENST').any():
            gene_names = pd.read_csv(scaff_path+'Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt', sep='\t', names=['ProbeID', 'G', 'Name'])
            del gene_names['G']
            df['ProbeID'] = df['ProbeID'].str.split('.', expand=True)[0]
            df = df.merge(gene_names, how='left', on='ProbeID')

        elif df['ProbeID'].str.contains('ENSMUST').any():
            gene_names = pd.read_csv(scaff_path+'Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt', sep='\t', names=['ProbeID', 'G', 'Name'])
            del gene_names['G']
            df['ProbeID'] = df['ProbeID'].str.split('.', expand=True)[0]
            df = df.merge(gene_names, how='left', on='ProbeID')

        elif df['ProbeID'].str.contains('ENSMUSG').any():
            gene_names = pd.read_csv(scaff_path+'Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt', sep='\t', names=['T', 'ProbeID', 'Name'])
            del gene_names['T']
            df['ProbeID'] = df['ProbeID'].str.split('.', expand=True)[0]
            df = df.merge(gene_names, how='left', on='ProbeID')

        else:
            df['Name'] = df['ProbeID']

        #df = df.merge(gene_names, how='left', on='ProbeID') Move Name column to correct spot
        names = df['Name']
        del df['Name']
        df.insert(1, 'Name', names)
        
        df = df.replace(np.NaN, 0)
    
        df.to_csv(path_dir+str(accessionID)+"-counts.txt", sep='\t', index=False)

    else:
        df = pd.read_csv(path_dir+str(accessionID)+"-counts.txt", sep='\t')
        

    expr = df.copy(deep=True)
    expr = expr.drop(['ProbeID', 'Name'], axis=1)

    print("Normalizing")
    adata = AnnData(expr.T)
    if norm == 'cpm':
        sc.pp.normalize_total(adata, target_sum=1e6)
    if takeLog:
        sc.pp.log1p(adata, base=2)

    norm_df = pd.DataFrame(adata.X)
    norm_df = norm_df.T
    
    norm_df.insert(0, 'ProbeID', list(df['ProbeID']))
    norm_df.insert(1, 'Name', list(df['Name']))
    norm_df.columns = list(df.columns)        

    norm_df.to_csv(path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name), header=True, index=False,sep='\t')
    
    print("Done with rna_analyze")


""" make_idx: Generates idx.txt

"""

def make_idx(gpl):
    print('Starting make_idx')
    expr = path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name)

    ptr = []
    ids = []
    name = []
    desc = []
    pos = 0

    with open(expr, 'rb') as f:
        for line in f:
            if pos == 0:
                pos += len(line)
            else:
                ptr.append(pos)
                pos += len(line)
                split = line.decode("utf-8").split('\t')
                ids.append(split[0])
                name.append(split[1].split(':')[0])
                desc.append(':'.join(split[1].split(':')[1:]))
        f.close()

    with open(path_dir+str(accessionID)+'-%s-idx.txt'%gpl.name, 'w') as f:
        f.write('ProbeID\tPtr\tName\tDescription\n')
        for i in range(len(ids)):
            f.write('{}\t{}\t{}\t{}\n'.format(ids[i], ptr[i], name[i], desc[i]))
        f.close()
    print("Done with make_idx")

""" survival_ih: generate survival and ih files

"""
def survival_ih(gpl):
    meta_df = None
    ih_dict = {'ArrayID':[], 'ArrayHeader':[], 'ClinicalHeader':[]}
    for name, gsm in gse.gsms.items():
        gsm_meta = {}
        gsm_meta['ArrayId'] = name
        gsm_meta['time'] = ''
        gsm_meta['status'] = ''

        ih_dict['ArrayID'].append(name)
        ih_dict['ArrayHeader'].append(name)
        ih_dict['ClinicalHeader'].append(gsm.metadata['title'][0])

        for col, val in gsm.metadata.items():
            # Metadata columns
            if col == 'title' or '_ch' in col:
                # Only keeping cols that don't have protocol and label -> not needed
                if 'protocol' not in col and 'label' not in col and 'channel_count' not in col:
                    # Sometimes all metadata in characteristics column, so need to split (ex age: 17, gender: male)
                    if ':' in val[0]:
                        for v in val:
                            # Split on ':' so everything before ':' is a new col and after ':' is the value
                            new_val = v.split(': ')
                            gsm_meta['c '+ new_val[0]+'_ch1'] = new_val[1]
                    else:
                        gsm_meta['c '+ col] = val

        # If need to add extra survival information
        # Ex:
        # meta_df['c description'] = meta_df['description']

        if meta_df is None:
            meta_df = pd.DataFrame(gsm_meta)
        else:
            temp_df = pd.DataFrame(gsm_meta)
            meta_df = meta_df.append(temp_df, sort=False)

    print(gsm.metadata.keys())
    meta_df.to_csv(path_dir+str(accessionID)+'-%s-survival.txt'%gpl.name, index=False, sep='\t')

    ih_df = pd.DataFrame(ih_dict)
    ih_df.to_csv(path_dir + str(accessionID) + '-%s-ih.txt' % (gpl.name), index=False, sep='\t')

    print('Done with survival_ih')

""" main(): runs the program

Args:
    accessionID: GEO accession ID
    takeLog: True if need to take log of data

"""

for _, gpl in gse.gpls.items():
    
    survival_ih(gpl)
    
    if all_files:
        if experimentType == 'm':
            m_analyze(gpl)
        elif make_expr:
            rna_analyze(gpl)
        make_idx(gpl)
    else:
        quit(1)
