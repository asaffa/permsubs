print("""
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# author. github/asaffa
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Permutation analysis of methylation data
# 1. Permutation procedure
#
# description: randomly permutes phenotype labels of subjects, performs t-tests, 
#              stores absolute t vals for all probes in each permutation
#
#input: matrix of DNA methylation M-values
#output: HDF5 file containing the results
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
""")
import sys
import os
import subprocess
import argparse
import gc
import gzip
import time
import math
import random
import re
import pdb
import itertools as itperm
import multiprocessing as multipro
import logging
import numpy as np
import pandas as pd
from scipy import stats as sp
import tables

def load_data(s,c):
    file_name = s 
    global in_data
    try:
        in_data = open(file_name)
    except IOError as e:
        logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
    
    pd.options.display.max_columns = 101
    df = pd.read_csv(in_data, sep = ',', header = 0, index_col=0, chunksize=c,iterator=False)

    return df

#function for randomly permuting case-control labels            
def rand_perm(i, df, rel):
    nperms = i
    n = len(df.columns)    
    
    #add row for phenotypic/case-control status
    case_control = np.empty([(n)], dtype=int)
    midpnt = int(np.ceil((len(case_control)/2)))
    case_control[0:midpnt] = 1
    case_control[(midpnt):len(case_control)] = 0
    df.loc["case-control",:] = case_control
    logger.info(df.loc["case-control",:]) 
    logger.info("Dataframe dimensions: %i rows by %i columns" %(
    len(df.index), len(df.columns)))
            
    #unpaired -> randomly assign labels by individual
    if rel == False:
        labels = list(df.loc["case-control", :])
        all_perms = []
        for i in range(0, int(nperms)):
            random.shuffle(labels)
            all_perms.append([bool(item) for item in labels])
    
    #output
    unique = [list(x) for x in 
              set(tuple(x) for x in all_perms)]
    logger.info("%i random permutations performed, of which %i are unique" %(
    int(nperms), len(unique)))
    return all_perms

#perform for every probe in every permutation
def mp_ttest_probes_in_perm(t):
    start = time.time()
    len_df = len(mp_meth_matrix)
    pval = list() 
    
    for perm in t:
        tup_p = tuple(perm)
        #ttest diff methylation for each probe, return vector of p-vals
        for m in range(0,(len_df)):
            bool_sel = np.array(perm)
            controls = mp_meth_matrix[m, bool_sel]
            cases = mp_meth_matrix[m, ~bool_sel]
            #perform ttest
            res = sp.ttest_ind(cases,controls,0,nan_policy="omit")
            #add to complete list of p-vals
            pval.append(abs(res[0]))
   
    logger.info('%s: finished' %(multipro.current_process().name))
    #print time-taken 
    elapsed = time.time() - start
    m,s = divmod(elapsed, 60)
    h,m = divmod(m, 60)
    logger.info("mp_ttest_probes_in_perm_unp: time taken: %d:%02d:%02d" %(h, m, s))
    logger.warn('Free memory after worker finished: {m} MB'.format(
                        m = free_memory()))
    #return complete list of p-vals (for all perms and all probes for this worker)
    return pval

def mp_write_results(r, c, n):
    #get individual results back from pooled results
    start = time.time()
    
    #num probes
    len_r = len(r)
    #init results table        
    names = np.array(r[0:(len_r)])
    global pvals_tab
    pvals_tab = pd.DataFrame(np.zeros((len_r, c),dtype='float32'))
    pvals_tab.index = names
    
    #counter for perms
    count_p = -1  
    #for each set of results from each worker
    for result in p_results:
        #for each of the lists of results 
        #(i.e. multiple permutation tests all flattened together)
        for j in range(0, len(result)):
            pval_temp = pd.DataFrame(np.zeros((len_r,np.int(c / n)),dtype='float32'))
            pval_temp.index = names
            #separate into individual perms and store as cols in results table
            for i in range(0,len(result[j]),(len_r)):
                count_p += 1
                #store row
                temp_row = result[j][i:(i + len_r)]
                pval_temp.ix[:,count_p] = temp_row
                
            #write to results table
            start_col = j * (c/n)
            end_col = j * (c/n) + (c/n)
            pvals_tab.ix[:, start_col : end_col] = pval_temp
            
    #print time-taken 
    elapsed = time.time() - start
    m,s = divmod(elapsed,60)
    h,m = divmod(m,60)
    logger.info("WriteResults: time taken: %d:%02d:%02d \n" %(h,m,s))  

# set up worker pool
def mp_worker_pool(d, n):
    len_d = len(d)
    quo_task_div, rem_task_div = divmod(len_d, n)
    task_list = list()
    start_it = 0
    
    for i in range(0, n):
        if i < rem_task_div:
            inc_it = quo_task_div
        else:
            inc_it = (quo_task_div - 1)
        task_list.append(d[start_it:(start_it + inc_it + 1)])
        start_it += inc_it + 1
    
    return task_list

def mp_main(ts, n):
    start = time.time()
    #number of processors
    num_proc = n
    mp_pool = multipro.Pool(num_proc)
    
    p_results = list()
    pool_arr = mp_pool.map(mp_ttest_probes_in_perm, ts)
    p_results.append(pool_arr)
    
    mp_pool.close()
    mp_pool.join()
    
    #print time-taken 
    elapsed = time.time() - start
    m,s = divmod(elapsed, 60)
    h,m = divmod(m, 60)
    logger.info("mp_main: time taken: %d:%02d:%02d \n" %(h, m, s))
    return p_results

def main():
    #these variables global so that each subprocess can share same memory
    global mp_meth_df, mp_meth_matrix, all_perms, ttest, cols, p_results, pvals_tab, store, logger
    logger = multipro.log_to_stderr(logging.INFO)
    logger.warn('Initial free memory: {m} MB'.format(m = free_memory()))
    
    #parse arguments from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("numcores", type = int)
    parser.add_argument("numperms", type = int)
    parser.add_argument("batchperms", type = int)
    parser.add_argument("chunksz", type = int)
    parser.add_argument("inputcsv")
    parser.add_argument("outputh5")
    args = parser.parse_args()

    start = time.time()
    #load methylation data
    mp_meth_df = load_data(args.inputcsv,2) 
    mp_meth_df = mp_meth_df.get_chunk()
    logger.info("First row of data (sorted by column/sample name) : ")
    logger.info(mp_meth_df.head(1))
    
    #generate permutations
    all_perms = rand_perm(args.numperms,mp_meth_df,False)
    num_proc = args.numcores
    num_perms = len(all_perms)
    
    #create hdf data store
    store = pd.HDFStore(args.outputh5)
    
    #########
    # split #
    #########
    #number of permutation to run at a time
    inc = args.batchperms
    if (num_perms % inc != 0) | (inc % num_proc != 0):
        raise Exception("batch size must divide both permutations and num procs exactly")
    logger.info("permutations batch size: {c}".format(c = inc))
    
    mp_meth_df = load_data(args.inputcsv,args.chunksz)
    
    for chunk in mp_meth_df:
        row_names = chunk.index.tolist()
        mp_meth_matrix = chunk.iloc[:,:].values
        logger.info("Dataframe dimensions: %i rows by %i columns" %(
        len(chunk.index),len(chunk.columns)))
        logger.warn('Free memory after loading data: {m} MB'.format(m = free_memory()))    
        for n in range(0, num_perms, inc):
            print("perms: %s" % (n))
            #logger.info("chunk = {m}".format(m=int(
            #                        (n / float(num_perms)) * int(num_perms / float(inc)))))
            #split between cores
            task_list = mp_worker_pool(all_perms[n:(n + inc)], num_proc)
            #########
            # apply #
            #########
            p_results = mp_main(task_list, num_proc)
            ###########
            # combine #
            ###########
            logger.warn('Free memory after all workers done: {m} MB'.format(
                                    m = free_memory())) 
            mp_write_results(row_names, inc, num_proc)
            hdf_node = str(round((n / float(num_perms)) * (num_perms / float(inc))))
            logger.warn('Free memory after combining results from workers: {m} MB'.format(
                                    m = free_memory())) 
            print("hd5 node: %s" %(hdf_node))
            store.append(hdf_node, pvals_tab, min_itemsize={ 'index' : 20 })
            logger.warn('Free memory after writing results to h5 file: {m} MB'.format(
                                    m = free_memory())) 
            gc.collect()
    
    #close hdf data store
    store.close() 
    
    #print time-taken 
    elapsed = time.time() - start
    m,s = divmod(elapsed,60)
    h,m = divmod(m,60)
    logger.info("main: time taken: %d:%02d:%02d" %(h,m,s))

def free_memory():
        total = 0
        a = list()
        p = subprocess.Popen('free -m', shell=True, stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
                a.append(line)
        total = int([x for x in a[1].split()][6])
        return total

def free_memory_mac():
        total = 0
        a = list()
        p = subprocess.Popen('/usr/bin/vm_stat', shell=True, stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
                a.append(line)
        tempint = int(a[1].strip(' \t\n\rPagesfree:.'))
        total = ((tempint * 4096) / 1024) / 1024
        return total

if __name__ == "__main__":
        main()



