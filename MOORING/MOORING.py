import numpy as np


def eoft(trmat):
    
    try:
        Trmat = trmat.data
        Trmat[trmat.mask] = np.nan
        trmat = Trmat
    except:
        None
    
    trmat = trmat.T-np.nanmean(trmat,axis=1)
    trmat = trmat.T
    
    mcov = np.dot(trmat,trmat.T)
    
    evals,evecs = np.linalg.eig(mcov)
    
    evecs_norm = evecs/np.linalg.norm(evecs,axis=0)
    
    evals_perc = evals/np.sum(evals)
    
    amp = np.dot(evecs_norm.T,trmat)
    
    return evals_perc,evecs_norm,amp
    