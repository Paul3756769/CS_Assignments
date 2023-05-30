# This file contains all functions i need for data analysis
import numpy as np
import matplotlib.pyplot as plt
from numba import njit


# Jackknife analysis

#@njit
def jackknife(data_raw, func,args, Blockcount = 0, t_int=0.5):
    # plausibility check for func and args
    if  callable(func)==False:
        raise ValueError("func is not callable")
    if isinstance(args, tuple)==False:
        raise ValueError("args is not a tuple")
    if t_int < 0:
        #raise ValueError("t_int is negative")
        t_int = 0.5
    step = int(2*t_int)+1

    data = data_raw[::step]

    blockcount = Blockcount
    n_raw = len(data)
    
    if Blockcount == 0:
        blockcount = n_raw
   
    blocklength = int(n_raw/blockcount)

    # truncate data
    n = blockcount*blocklength
    data = data[(n_raw-n):n_raw]

    value = func(data,*args)

    jk = np.zeros(blockcount)
    for i in range(blockcount):
        d_jk = np.delete(data, range(i*blocklength, (i+1)*blocklength))
        jk[i] = func(d_jk,*args)
    #    print(np.mean(data[i*blocklength:(i+1)*blocklength]),np.mean(d_jk))
        #print(jk[i],value,len(d_jk),len(data))	

    std_err = np.sqrt((blockcount-1)/blockcount*np.sum((jk-value)**2))
    return value, std_err

@njit
def err_t_int(t_int,k_max,N):
    return t_int*np.sqrt(2*(2*k_max+1)/N)

@njit
def crit(t,k):
    return (k >= 6*t)


@njit
def Autocorrelation_series(ts):
    N = len(ts)

    A = np.zeros(N)
    t = np.zeros(N)

    t_int = 0.5

    k_crit = N-1
    k_max = N

    mean_ts = np.mean(ts)
    var_ts = np.var(ts)

    check = True

    for k in range(1,N):
        if k is k_max:
            break
        ts1 = ts[k:]
        ts2 = ts[:-k]
        a1 = np.mean(ts1*ts2)
        a2 = mean_ts**2
        a = (a1-a2)/var_ts
        A[k] = a
        t_int += a
        t[k] = t_int

        if t_int<0 and check:
            k_crit = 0
        if crit(t_int,k) and check:
            k_crit = k
            k_max = get_kmax(k_crit,N)
            check = False

        
    A_reshaped = A[:k_max]
    t_reshaped = t[:k_max]


    return A_reshaped, t_reshaped, k_crit
@njit
def get_kmax(k_crit,N):
    return min(5*k_crit,N)

# binning analysis

def autocorrelation_time_binning(TS):
    base = 2
    
    len_TS = len(TS)
    ex0 = 7
    while len_TS < base**ex0:
        ex0 -=1

    ex = ex0
    
    while len_TS > base**(ex+1):
        ex += 1
    pow2 = base**ex
    quotient = len_TS//pow2
    remainder = len_TS%pow2
    ts = TS[remainder:]

    binsizes = np.array([quotient*base**i for i in range(max(0,ex0-4),max(0,ex-4))])

    variance_corr = np.var(ts)
    variances_uncorr = np.zeros(len(binsizes))


    for i in range(len(binsizes)):
        binsize = binsizes[i]
        bin_count = int(len(ts)/binsize)
        bin_means = np.zeros(bin_count)
        for j in range(bin_count):
            bin = ts[j*binsize:(j+1)*binsize]
            bin_means[j] = np.mean(bin)
        variances_uncorr[i] = np.var(bin_means)
        
    



    N = len(ts)
    
    #t_ints = variances_uncorr/variance_corr/2*binsizes
    
    t_ints = np.zeros(len(binsizes))
    for i in range(len(binsizes)):
       try: 
        t_ints[i] = variances_uncorr[i]/variance_corr/2*binsizes[i]
       except OverflowError:
        t_ints[i] = 0    
    bincounts = N/binsizes
    
    return bincounts,t_ints,binsizes
    





# bivariate ts

@njit
def rnd_gauss_0_1():
    return np.random.normal(0,1)
@njit
def bivariate_ts(p,length):

    ts = np.zeros(length)
    ts[0] = rnd_gauss_0_1()
    for i in range(1,length):
        r = rnd_gauss_0_1()
        ts[i] = ts[i-1]*p + r*np.sqrt(1-p**2)
    t_int = 0.5 * (1+p)/(1-p)
    return ts, t_int


