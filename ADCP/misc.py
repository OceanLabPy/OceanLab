import numpy as np
import scipy.interpolate as scint

def nans(size):
    return np.zeros(size)*np.nan
    

def adcp_binning(ADCP,latadcp,lonadcp,latctd,lonctd):
    '''
    Important to use on geostrophic velocity reference,
    this function calculates the bin average velocity,
    profile from ADCP data between the CTD stations.
    To interpolates ADCP data into regular spaced depths,
    of 1 meter use adcp_equalpress().
    
    Usage:
         adcp_binning(2D-np.array,1D-np.array,1D-np.array,
         1D-np.array,1D-np.array) -> 2D-np.array
    '''
    latadcp,lonadcp = np.array(latadcp),np.array(lonadcp)
    latctd,lonctd   = np.array(latctd),np.array(lonctd)

    args = argdistnear(latctd,lonctd,latadcp,lonadcp)
    
    ADCP_m = []
    for inf,sup in zip(args[:-1],args[1:]):
        ADCP_m.append(np.mean(ADCP[:,inf:sup],axis=1))
    ADCP_m = np.vstack(ADCP_m).T
    
    return ADCP_m


def adcp_equalpress(ADCP,PRESS,kind='linear'):
    pp = np.arange(0,PRESS.max()+1)
    ADCP_pp = []
    for col in np.arange(ADCP.shape[1]):
         f = scint.interp1d(PRESS[:,col],ADCP[:,col],kind=kind)
         ADCP_pp.append(f(pp[pp>=PRESS.min()]))
    ADCP_pp = np.vstack(ADCP_pp).T
    
    if PRESS.min()>0:
        ADCP_pp = np.vstack([nans((PRESS.min(),ADCP.shape[1])),
                            ADCP_pp])
    PP = np.tile(pp,(ADCP.shape[1],1)).T
    return ADCP_pp,PP
    
    