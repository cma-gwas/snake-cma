import pandas as pd
import numpy as np
from math import sqrt
from scipy.stats import chi2
import argparse
import os
import time
import datetime

def initial_logging(log, inputs, inputs_Null, output, args_split, args_meta_summary, args_meta_no_correction, args_out):
    assert args_split is not None, "Split Index is required.."
    assert args_meta_summary or args_meta_no_correction, "Choose a flag for meta-analysis procedure. Currently avaiable options are: --meta_summary and --meta_no_correction."
    assert args_out is not None, "Output path is required.."
    with open(output, "w") as f:
        f.write(log % (datetime.datetime.now().strftime("%A - %d %B %Y - %H:%M:%S"), str(os.getcwd())))
        for i in inputs:
            if inputs[i]!= None and type(inputs[i])!=bool:
                f.write("\t--%s %s\n" % (str(i), str(inputs[i])))
            elif inputs[i]!=inputs_Null[i] and type(inputs[i])==bool:
                f.write("\t--%s \n" % str(i))

def import_results(path_prefix_, split_index_, cols_=['A1FREQ','BETA', 'SE'], cols__=['CHROM','GENPOS','ID','ALLELE0','ALLELE1']):
    # split_index_ > 1
    df_temp = pd.read_csv(path_prefix_ % 1, usecols= cols__ + cols_, sep=" ", header=0)
    beta_ = np.array(df_temp['BETA'], dtype=np.float32).reshape((-1,1,1))
    se_ = np.array(df_temp['SE'], dtype=np.float32).reshape((-1,1,1))
    n_ = 2 * len(df_temp)
    f_ = np.array(df_temp['A1FREQ'], dtype=np.float32) * n_
    for i in range(2,split_index_+1):
        df_temp_2 = pd.read_csv(path_prefix_ % i, usecols= cols_, sep=" ", header=0)
        beta_ = np.hstack((beta_,np.array(df_temp_2['BETA'], dtype=np.float32).reshape((-1,1,1))))
        se_ = np.hstack((se_,np.array(df_temp_2['SE'], dtype=np.float32).reshape((-1,1,1))))
        f_ = f_ + np.array(df_temp_2['A1FREQ'],dtype=np.float32) * 2 * len(df_temp_2)
        n_ = n_ + 2 * len(df_temp_2)
    df_temp['A1FREQ'] = f_/n_
    df_temp['A1FREQ'] = np.round(np.array(df_temp['A1FREQ']),6)
    return (beta_,se_,df_temp[cols__+['A1FREQ']])

def compute_cor_summary(beta_,se_,correction_bool=True):
    z_ = beta_ / se_
    alpha_ = np.array([[(z_[:,i,0] * z_[:,j,0])/np.abs(z_[:,i,0] * z_[:,j,0]) for i in range(len(z_[0]))] for j in range(len(z_[0]))]).transpose(2,0,1)
    if correction_bool:
        return alpha_ * np.repeat(np.array([[round(np.corrcoef(z_[:,i,0],z_[:,j,0])[0,1],13) for i in range(len(z_[0]))] for j in range(len(z_[0]))], dtype=np.float32)[np.newaxis, :, :], len(se_), axis=0)
    else:
        return np.repeat(np.array([[round(np.corrcoef(z_[:,i,0],z_[:,j,0])[0,1],13) for i in range(len(z_[0]))] for j in range(len(z_[0]))], dtype=np.float32)[np.newaxis, :, :], len(se_), axis=0)

def compute_cov(corr_matrix, se_):
    return se_ * corr_matrix * se_.transpose(0,2,1) # returns M x K x K numpy array

def compute_weights(cov_matrix, se_, split_index_):
    # M = len(se_) is the number of Markers
    # corr_matrix is of dimension M x K x K where K = split_index
    # se_ and beta_ are matrices of dimension M x K x 1 where K = split_index

    e_temp_1 =  np.repeat(np.ones(split_index_, dtype=np.float32).reshape((split_index_,1))[np.newaxis, :, :], len(se_), axis=0)
    e_temp_2 =  np.repeat(np.ones(split_index_, dtype=np.float32).reshape((1,split_index_))[np.newaxis, :, :], len(se_), axis=0)
    M = len(se_)

    ## We compute M many inverse of Cov Matrices (K x K). Because these are independent, we have split this step into two for resource efficiency purposes..
    inv_cov_matrix_1 = np.linalg.inv(cov_matrix[:int(len(se_)/2)]) # returns M x K x K numpy array
    inv_cov_matrix_2 = np.linalg.inv(cov_matrix[int(len(se_)/2):]) # returns M x K x K numpy array
    inv_cov_matrix= np.concatenate([inv_cov_matrix_1,inv_cov_matrix_2], axis=0)
    return np.matmul(e_temp_2, inv_cov_matrix).reshape(-1,split_index_,1)/np.matmul(np.matmul(e_temp_2, inv_cov_matrix), e_temp_1)

def compute_stats(weights, cov_matrix, beta_, se_, split_index_, df_):
    beta_final_ = np.matmul(beta_.reshape(-1,1,split_index_),weights).reshape(-1,1)
    se_final_ = np.sqrt(np.matmul(np.matmul(weights.transpose(0,2,1),cov_matrix),weights).reshape(-1,1))
    df_['BETA'] = np.round(beta_final_,6)
    df_['SE'] = np.round(se_final_,6)
    df_['CHISQ'] = (np.array(df_['BETA'])/np.array(df_['SE']))**2
    df_['LOG10P'] = -1 * np.log10(chi2.sf(df_['CHISQ'], 1, loc=0, scale=1))
    df_['P'] = 10**((-1) * np.array(df_['LOG10P']))
    df_ = df_.set_index('CHROM')
    df_['CHISQ'] = df_['CHISQ'].apply('{:,.3e}'.format)
    df_ = df_[['GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'CHISQ', 'LOG10P', 'P']]
    print('Done..')
    return df_

def meta_nocorrection(path_prefix_, split_index_):
    if split_index_ == 1:
        df_ = pd.read_csv(path_prefix_ % 1, usecols= ['CHROM','GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'CHISQ', 'LOG10P'], sep=" ", header=0)
        df_ = df_.set_index('CHROM')
        df_['P'] = -1 * np.log10(chi2.sf(df_['CHISQ'], 1, loc=0, scale=1))
        return df_
    else:
        beta_,se_,df_ = import_results(path_prefix_, split_index_)
        print('Results are imported..')
        M = sum(1 for line in open(path_prefix_ % 1)) - 1
        corr_matrix_ = np.repeat(np.identity(split_index_)[np.newaxis,:,:],M,axis=0)
        print('Correlation Matrix is computed..')
        cov_matrix_ = compute_cov(corr_matrix_, se_)
        print('Covariance Matrix is computed..')
        del corr_matrix_, M
        weights_ = compute_weights(cov_matrix_, se_, split_index_)
        print('Weights computed..')
        return compute_stats(weights_, cov_matrix_, beta_, se_, split_index_, df_)

def meta_summary(path_prefix_, split_index_):
    if split_index_ == 1:
        df_ = pd.read_csv(path_prefix_ % 1, usecols= ['CHROM','GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE', 'CHISQ', 'LOG10P'], sep=" ", header=0)
        df_ = df_.set_index('CHROM')
        df_['P'] = -1 * np.log10(chi2.sf(df_['CHISQ'], 1, loc=0, scale=1))
        return df_
    else:
        beta_,se_,df_ = import_results(path_prefix_, split_index_)
        print('Results are imported..')
        corr_matrix_ = compute_cor_summary(beta_,se_)
        print('Correlation Matrix is computed..')
        cov_matrix_ = compute_cov(corr_matrix_, se_)
        print('Covariance Matrix is computed..')
        del corr_matrix_
        weights_ = compute_weights(cov_matrix_, se_, split_index_)
        print('Weights computed..')
        return compute_stats(weights_, cov_matrix_, beta_, se_, split_index_, df_)


