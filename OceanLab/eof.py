import numpy as np

# functions
#=========================================
# COMPUTE EOFS
#=========================================
def eoft(trmat,nm=None):
    '''
    evals_perc,evecs_norm,amp = eoft(trmat)

    Computes eof in time (or depth) for general cases
    normally trmat is a matrix containing each time series
    (or vertical profiles) as a row.
    That is, for N stations having m data points,
    trmat is a N by m matrix.
    '''

    #if is masked the data
    try:
        Trmat = trmat.data
        Trmat[trmat.mask] = np.nan
        trmat = Trmat.copy()
    except:
        None

    #% demeans trmat prior to doing eof
    trmat = (trmat.T-np.nanmean(trmat,axis=1)).T
    #% computes zero-lag cross-covariance matrix
    mcov = np.dot(trmat,trmat.T)/(trmat.shape[1]-1)
    #% computes eigenvalues and eigenvectors
    evals,evecs = np.linalg.eig(mcov)
    #find the descending order of eigenvalues
    ind = np.argsort(evals)[::-1]
    #% normalize eigenvectors
    evecs_norm = evecs[:,ind]/np.linalg.norm(evecs[:,ind],axis=0)
    #% sort eigenvalues and computes percent variance explained by each mode
    evals_perc = evals[ind]/np.sum(evals)
    #% computes amplitude functions
    amp = np.dot(evecs_norm.T,trmat)

    #if was chosen a number of eigenvectors
    if nm!=None:
        evecs_norm = evecs_norm[:,:nm]
        evals_perc = evals_perc[:nm]
        amp        = amp[:nm,:]

    return evals_perc,evecs_norm,amp

#=========================================
# FILL GAPS WITH EOFS
#=========================================
def my_eof_interp(M,nmodes,errmin=1e-15,repmax=None):

    """
    vi = my_eof_interp(M,nmodes)

    This function try to fill gappy data using its EOFs

    INPUT:
    -> M: the data Matrix (generaly the diferent stations are on the columns);
    -> nmodes: the number of statiscally significant EOF modes (to be used for the series reconstruction)

    The method -- based on Beckers & Rixen (2003, JAOTech)

    The algorithm proposed by the authors is as follows.  (see Eqs. (13) and (14) from Beckers & Rixen)
        (I) Put zero in the GAPPY positions ("unbiased INITIAL GUESS, by assuming that the removed mean was unbiased).
        (II) Then you can have a first estimative of the EOFs.
        (III) You can reconstruct the series using the statiscally significant modes (e.g.: determined using Monte Carlo simulations as in Preisendorfer [1988] -- OBS: the authors propose another method using cross validation during the iterative process)
        (IV) Substitute the reconstructed series ONLY in the GAPPY positions in your original matrix.
        (V) Run steps (II), (III) and (IV) iteractively till its convergence (i.e. compare the reconstructed series at time t+1 with that at time t).

    %%%%% OBS: my_eof_interp doesn't work for filling gaps simultaneously in all depths.
    """

    Mmean = np.nanmean(M,axis=1)

    M = M.T - np.nanmean(M,axis=1)
    M = M.T

    f = np.isnan(M)
    M[f] = 0 #%% initial guess in gappy positions

    rep = 0
    while True:
       print(rep+1)


       Ud,Dd,Vd = np.linalg.svd(M,'econ')
       Dd = np.diagflat(Dd)

       #% compute the explained variance (total variance is the sum of the Dd elements squared -total energy-)
       evd = np.diag(Dd)**2 #%% remember that the singular values (Dd) are the sqrt of
       evd = Dd/np.sum(Dd)  #%%    the eigenvalue of the cov. matrix.


       #%% truncated EOF reconstruction
       Ud = Ud[:,0:nmodes]
       Vd = Vd[:,0:nmodes]
       Dd = Dd[0:nmodes,0:nmodes]
       #Dd = Dd[0:nmodes]

       #% reconstructing the series -- truncated reconstruction (SVD: Mrec = Ud * Dd *Vd')
       Xa = np.dot(np.dot(Ud,Dd),Vd.T)

       Mold = M.copy() #% saving for compute the error

       M[f] = Xa[f]

       #% testing the error
       if rep>=1:
          err = np.nanmean(np.abs(M - Mold))/(np.nanmax(M)-np.nanmin(M))
          print(err)
          if err < errmin:
             break
       if repmax!=None:
            if rep>=repmax:
                break
       rep += 1

    #%% Saving the last the series with the EOF reconstructed series in the gappy positions;
    vi = (M.T+Mmean).T

    return vi

#=========================================
# PERFORM COMPLEX EOF
#=========================================
from scipy.signal import hilbert
import scipy.linalg as la
import numpy as np
from dask import delayed
from scipy.signal import hilbert
import scipy.linalg as la
import numpy as np

def ceof(data, nkp = 10):
    ''' Complex (Hilbert) EOF
    First written in MATLAB and found in Prof. Daniel J. Vimont webpage 
    (https://www.aos.wisc.edu/~dvimont/matlab/Stat_Tools/complex_eof.html)
    
    Translated to Python by Felipe Vilela da Silva @ IMAS, UTas, AU on 15/Mar/2021
    ==============================================================================
    Inputs:
    -> data: original data set [time, space]
    -> nkp:  number of modes to output (default = 10)

    ==============================================================================
    Outputs:
    -> lamda: eigenvalues
    -> load:  First nkp Complex Loadings or eigenvectors [space, nkp]
    -> pcs:   First nkp Complex Principal Components or amplitudes [time, nkp]
    -> per:   percent variance explained (real eigenvalues)
    ==============================================================================
    Version 2.0.0 by Felipe Vilela da Silva on 25/May/2021. 
        Now, it is possible to input data with NaN values. 
        The quantity of lambda is related the amount of not Nan values in data
    '''
    load_real = np.zeros([data.shape[1], nkp])*np.nan
    load_imag = np.zeros([data.shape[1], nkp])*np.nan
    # It is necessary to remove the nan values of the matrix to solve the eigenvalue problem
    nan_values = np.isnan(data[0,:]) # We can just look at each coordinate along a single time
    data = data[:,~nan_values]   # Then, we remove all these coordinates in all of the occurences
    
    ntim, npt = data.shape
    
    # Hilbert transform: input sequence x and returns a complex result of the same length
    print('1: Performing Hilbert transform')
    data_hilbert = hilbert(data)
    # Compute the covariance matrix in the Hilbert transform
    print('2: Computing covariance matrix')
    c = delayed(np.dot)(data_hilbert.conjugate().T, data_hilbert).compute()/ntim
    print('3: Solving the eigenvalue problem')
    lamda, loadings = delayed(la.eig)(c).compute() # lamda: eigenvalue, loadings: eigenvectors
    
    l = lamda.conjugate().T; k = np.argsort(l)
    lamda, loadings = np.flip(l[k]), np.fliplr(loadings[:,k])
    loadings = loadings[:,:nkp]
    
    per = lamda.real*100/np.sum(lamda.real)
    pcs = np.dot(data_hilbert,loadings)
    
    load_real[~nan_values,:] = loadings.real.copy()
    load_imag[~nan_values,:] = loadings.imag.copy()
    load = load_real + 1j*load_imag
    print('Done! \U0001F600')
    
    return lamda, load, pcs, per

def amplitude_phase(evecs, amp):
    ''' Complex (Hilbert) EOF
    First written in MATLAB and found in the webpage below 
    (https://www.jsg.utexas.edu/fu/files/GEO391-W11-CEOF.pdf)
    
    Translated to Python by Felipe Vilela da Silva @ IMAS, UTas, AU on 15/Mar/2021 
    ===========================================================================
    Inputs:
    -> evecs: First nkp Complex Loadings or eigenvectors [space, nkp]
    -> amp:   First nkp Complex Principal Components or amplitudes [time, nkp]

    ===========================================================================
    Outputs:
    -> sp_amp:   Spatial amplitude [space, nkp]
    -> sp_phase: Spatial phase [space, nkp]
    -> t_amp:    Temporal amplitude [time, nkp]
    -> t_phase:  Temporal phase [time, nkp]
    ===========================================================================
    '''
    # Spatial amplitude
    sp_amp = pow(np.multiply(evecs, np.conj(evecs)),0.5)
    theta = np.arctan2(evecs.imag, evecs.real)
    # Spatial phase
    sp_phase = np.divide(np.multiply(theta, 180), np.pi)

    # Temporal amplitude
    t_amp = pow(np.multiply(amp, np.conj(amp)), 0.5)
    # Temporal phase
    phit = np.arctan2(amp.imag, amp.real)
    t_phase = np.divide(np.multiply(phit, 180), np.pi)
    
    return sp_amp, sp_phase, t_amp, t_phase 

def reconstruct_ceof(amp, modes, n, day):
    ''' Reconstrucion of CEOF modes individually following Majumder et al. (2019)
    
    Written by Felipe Vilela da Silva @ IMAS, UTas, AU on 15/Apr/2021 
    ===========================================================================
    Inputs:
    -> amp:    amplitude or coefficient of expansion [time, n]
    -> modes:  eigenvector or loading [lat, lon, n]
    -> n:      mode to be reconstructed. 0 is the first
    -> day:    day to be reconstructed. 0 is the first
    ===========================================================================
    Output:
    -> Rec_ceof: Reconstruction of the CEOF field [time, lat, lon]
    ===========================================================================
    '''   
    
    # Majumder et al (2019) compute the reconstructed CEOF field as the real part of the multiplication between
    # the coefficient of expansion (i.e., amplitude) and the complex conjugate of the loading (i.e., mode)
    Rec_ceof = amp[day,n]*np.conj(modes[:,:,n])
    return Rec_ceof.real

