import numpy as np

from OceanLab.utils import *

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
        ADCP_m.append(np.nanmean(ADCP[:,inf:sup],axis=1))
    ADCP_m = np.vstack(ADCP_m).T

    return ADCP_m


def mdr(geostrophic,observation,p_ini=100,p_fin=200):
        '''
        This function calibrate the normal geostrophic velocity
        section from observed normal velocity data. It must be in
        the same shape to work.
        p_ini and p_fin defines the index of range of reference levels.
        The method is based on find the level which once used to reference
        the geostrophic normal section data reduces to maximum the mean
        squared error.

        Usage:
            mdr(2D-array,2D-array,integer,integer) --> 2D-array


        Dependence:
            rmse()

        '''

        p_ini = int(p_ini)
        p_fin = int(p_fin)

        vsec = geostrophic.copy()
        adcp = observation.copy()

        if p_fin >= adcp.shape[0]:
            print('The given depth %.i is out of the max data depth.'%(p_fin))
            p_fin = adcp.shape[0]-1
            print('Using now %.i'%(p_fin))


        #Set pref as an empty list
        pref = []
        #reads each column
        for col in np.arange(0,vsec.shape[1]):
                #set RSME as an empty list
            rsmes = []
            for lvl in np.arange(p_ini,p_fin):
                reff = adcp[lvl,:]-vsec[lvl,:]
                vsecref = vsec+reff
                rsmes.append(rsme(vsecref[0:adcp.shape[0],col],
                                     adcp[:,col]))

            if np.all(np.isnan(rsmes)):
                pref.append(p_fin+1)
            else:
                pref.append(np.nanargmin(rsmes)+p_ini)

        pref = np.array(pref).astype('int')
        cols = np.arange(0,adcp.shape[1])

        adcpvals = adcp[pref,cols]
        vsecvals = vsec[pref,cols]
        adcpvals[pref==p_fin+1] = np.nan
        vsecvals[pref==p_fin+1] = np.nan

        ref = adcpvals-vsecvals
        vmdr = vsec+ref

        return vmdr,pref
