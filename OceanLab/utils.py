import numpy as np
import scipy.signal as sg
import xarray as xr

import dask
from dask.distributed import Client, progress

##### User functions
#=============================================================================
# NEAREST DISTANCE
#=============================================================================
def argdistnear(x,y,xi,yi):
    '''
    This function finds the index to nearest points in (xi,yi) from (x,y)
    ======================================================================

    USAGE:
    x,y = [5,1,10],[2,6,3]
    xi,yi = np.linspace(0,19,20),np.linspace(-5,30,20)
    ind = argdistnear(x,y,xi,yi)

    INPUT:
       (x,y)   = points [list]
       (xi,yi) = series to search nearest point [list]

    OUTPUT:
       ind     = index of the nearest points
    ======================================================================
    '''
    
    idxs = [np.argmin(np.sqrt((xi-xx)**2 + (yi-yy)**2)) for xx,yy in zip(x,y)]
    idxs = np.array(idxs)
    return idxs
#=============================================================================

#=============================================================================
# LOW PASS FILTER
#=============================================================================
def meaneddy(prop,days=60,ndim=1,DataArray=False,timedim=None):

    """
    Apply a low-pass filter (scipy.signal.butter) to 'prop' and obtain the 
     mean and eddy components.
    ==========================================================================

    USAGE [1]:
    # np.array() one year, 17 lat x 13 lon domain
    Velocity = np.random.randn(365,17,13)
    Filtered, Residual = meaneddy(Velocity, days=10, ndim=3,
                                   DataArray=False,timedim=None)
    
    USAGE [2]:
    # xr.DataArray() one year, 17 lat x 13 lon domain
    Velocity = xr.DataArray(data=np.random.randn(365,17,13),
                             dims=["time","lat","lon"],
                             coords=dict(time=(["time"],range(0,365)),
                             lat=(["lat"],np.arange(-4,4.5,0.5)),
                             lon=(["lon"],np.arange(1,7.5,0.5))))
    Filtered, Residual = meaneddy(Velocity, days=10,
                                   DataArray=True,timedim=["time"])
    
    INPUT:
       prop      = 1, 2 or 3D array to be filtered
       days      = number of days to set up the filter
       ndim      = number of dimensions of the data
                    [only used for DataArray=False, max:3]
       DataArray = True if prop is in xr.DataArray format
       dim       = name of time dimension to filter
                    [only used for DataArray=True]
   
    OUTPUT:
       m_prop    = mean (filtered) part of the property
       p_prop    = prime part of the property, essentially prop - m_prop
    ==========================================================================
   """

    # creating filter
    def timefilter(prop,filtdays=60):
        filt_b,filt_a = sg.butter(4,1./filtdays)
        return sg.filtfilt(filt_b,filt_a,prop)

    if DataArray:
        m_prop = xr.apply_ufunc(
                timefilter, # first the function
                prop,# now arguments in the order expected by 'butter_filt'
                input_core_dims=[timedim],  # list with one entry per arg
                output_core_dims=[timedim],  # returned data
                kwargs={'filtdays':days},
                vectorize=True,  # loop over non-core dims
                dask='vectorized')
        
        p_prop = prop - m_prop
        
    elif ndim ==3:
        ti,lt,ln = prop.shape
        prop = prop.reshape(ti,lt*ln)

        m_prop = []

        for tot in prop.T:
        # filtered series (mean)
            m_prop.append(timefilter(tot,filtdays=days))

        m_prop = np.array(m_prop)
        m_prop = m_prop.T

        p_prop = prop - m_prop

        m_prop = m_prop.reshape(ti,lt,ln)
        p_prop = p_prop.reshape(ti,lt,ln)

    elif ndim ==2:

        ti,lt = prop.shape

        m_prop = []

        for tot in prop.T:
            # filtered series (mean)
            m_prop.append(timefilter(tot,filtdays=days))

        m_prop = np.array(m_prop)
        m_prop = m_prop.T

        p_prop = prop - m_prop

        m_prop = m_prop.reshape(ti,lt)
        p_prop = p_prop.reshape(ti,lt)

    elif ndim ==1:
        m_prop = timefilter(prop,filtdays=days)
        p_prop = prop - m_prop

    return m_prop,p_prop
#=============================================================================

##### Functions for relative imports
# =============================================================================
# KERNEL FOR PARALLEL COMPUTING
# =============================================================================
def _parallel_client(cpu_params=dict(tpw=2,nw=4,ml=7.5)):
    """
    Create client kernel for parallel computing
    ====================================================
    INPUT:
        -> cpu_params: dict containing floats with keys
            -> tpw: threads_per_worker
            -> nw: n_workers
            -> ml: memory_limit per worker [GB]
    OUTPUT:
        -> client: configuration of parallel computing
    ====================================================
    """

    client = Client(threads_per_worker=cpu_params['tpw'], 
                    n_workers=cpu_params['nw'], 
                    memory_limit=str(cpu_params['ml'])+'GB')
    return client
#=============================================================================



