import sys
import os.path
import gc
import argparse
import math
import itertools as itperm
import random
import re
import time
from timeit import Timer
import multiprocessing as multipro
import numpy as np
import pandas as pd
import tables
from scipy import stats as sp
########################
#load permutation data #
########################

def load_data(file_name):
    #global perms
    global h5
    try:
        in_data = open(file_name)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
    h5 = pd.HDFStore(file_name)
    #return perms
    return h5
    
def select_frame(n):
    global perms
    perms = h5.select(n)
    return perms

########################
# get p values         #
# for table of abs Ts  #
########################

def mp_sub_sample(t):
    global perms
    np.random.seed()
    idx_start = t[0]
    idx_fin = t[1]

    all_res = list()
    npop = len(perms)
    nreps = 1
    inc = int((1/float(100)) * npop)
    
    for p in range(idx_start, idx_fin):
        print "p: %s" %(p)
        pop = np.array(perms.ix[:,p])
        pop = sp.t.sf(pop,(ncase - 2)) * 2
        p_res = np.zeros(101)
        p_res[0] = 0
    
        for reps in xrange(0,nreps):
            np.random.shuffle(pop)
            min_p = np.zeros(101) * np.nan
            min_p[0] = 0
            #for each density, get the min p of subsample starting from 
            #position 0 in the pop up to current density 
            for dens in xrange(1,101):
                min_p[dens] = np.amin(pop[0:int((dens * inc) + 1)])
            #add results for this rep to the cumulative results
            p_res += min_p
        #divide cumulative results by reps to get mean values for each probe density    
        p_res = p_res/(float(nreps))
        all_res.append(p_res.tolist())
        del min_p
        del p_res
        del pop
    return all_res


# In[71]:

def mp_worker_pool(d, n):
    
    quo_task_div, rem_task_div = divmod(d, n)
    task_list = list()
    start_it = 0
    
    for i in xrange(0, n):
        if i < rem_task_div:
            inc_it = quo_task_div
        else:
            inc_it = (quo_task_div - 1)
        task_list.append([start_it,(start_it + inc_it + 1)])
        start_it += inc_it + 1
    print 'task list: %s' %(task_list)
    return task_list


# In[171]:

def mp_main(ts, n):
    start = time.time()
    num_proc = n
    mp_pool = multipro.Pool(num_proc)
    l_res = list()
    pool_arr = mp_pool.map(mp_sub_sample, ts)
    l_res.append(pool_arr)
    mp_pool.close()
    mp_pool.join()
    
    #print time-taken 
    elapsed = time.time() - start
    m,s = divmod(elapsed, 60)
    h,m = divmod(m, 60)
    print "mp_main: time taken: %d:%02d:%02d \n" %(h, m, s)    
    return l_res

def to_table(lr, n):
    
    ncol = 101
    dtypes = np.repeat('f8',ncol)
    names = ["{:01d}".format(x) for x in range(0,ncol)] 
    cols = zip(names,dtypes)
    tab = pd.DataFrame(np.zeros(n,dtype=cols))
    count = -1

    #get results back from each worker
    for worker in lr:
        for res in worker:
	    for perm in xrange(0,(len(res))):
                count += 1
                tab.iloc[count,:] = (res[perm])
    
    return tab

####################################
# write subsampling results to csv #
####################################

def write_file(file_name):
    global res_tab
    #append to file (or create new if none exists)   
    file_out = open(file_name, 'a')
    pd.DataFrame.to_csv(res_tab,file_out.name,sep = ',',mode='a',
    header = False, index = False)
    file_out.close()

def main():
    global h5, perms, num_proc, res_tab, ncase
    parser = argparse.ArgumentParser()
    parser.add_argument("numcores", type = int)
    parser.add_argument("inputcsv")
    parser.add_argument("ncase", type = int)
    parser.add_argument("st", type = int)
    parser.add_argument("fn", type = int)
    parser.add_argument("reps", type = int)
    parser.add_argument("outputcsv")
    args = parser.parse_args()
    ncase = args.ncase
    #start timer
    start = time.time()
    
    #hd5 file
    h5 = load_data(args.inputcsv)
    
    perms = select_frame("0")
    len_df = len(perms.columns)
    print "Dataframe dimensions: %i rows by %i columns" %(
            len(perms.index),len(perms.columns))
    task_list = mp_worker_pool(len_df, args.numcores)
    nreps = args.reps
    for reps in xrange(0,nreps):
        #p_res = pd.DataFrame(columns=xrange(0,101))
        p_res = pd.DataFrame(columns=[str(x) for x in xrange(0,101)])
        print "rep: %s" %(reps)
        np.random.seed()
        for df in xrange(args.st,args.fn):
            print "dataframe : %s" %(df)
            perms = select_frame(str(df))
            l_res = mp_main(task_list, args.numcores)
            p_res = p_res.append(to_table(l_res, len_df),ignore_index=True)
        #1/0
        res_tab = p_res.apply(lambda x: np.percentile(x,5),0)
        res_tab = pd.DataFrame(res_tab).T
        #1/0
        write_file(args.outputcsv)
    #print time-taken
    elapsed = time.time() - start
    m,s = divmod(elapsed, 60)
    h,m = divmod(m, 60)
    print "main: time taken: %d:%02d:%02d \n" %(h, m, s)    

if __name__ == "__main__":
    main()
