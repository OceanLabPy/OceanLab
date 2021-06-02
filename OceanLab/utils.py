import numpy as np
import scipy.signal as sg
import xarray as xr

##### Function
#=============================================================================
# NEAREST DISTANCE
#=============================================================================
def argdistnear(x,y,xi,yi):
    '''
    This function finds the index to nearest points in (xi,yi) from (x,y).
    
    usage:
    x,y = [5,1,10],[2,6,3]
    xi,yi = np.linspace(0,19,20),np.linspace(-5,30,20)
    ind = argdistnear(x,y,xi,yi)

    INPUT:
    --> (x,y): points [list]
    --> (xi,yi): series to search nearest point [list]
    
    Iury T.Simões-Sousa
    (IO-USP/ UMass-Dartmouth)
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
    Apply a low-pass filter (scipy.signal.butter) to 'prop' and obtain the mean and eddy components.
    
    usage [1]:
    Velocity = np.random.randn(365,17,13) # one year, 17 lat x 13 lon domain
    Filtered, Residual = meaneddy(Velocity, days=10, ndim=3, DataArray=False,timedim=None)
    
    usage [2]:
    Velocity = xr.DataArray(data=np.random.randn(365,17,13), dims=["time","lat","lon"],
coords=dict(time=(["time"],range(0,365)), lat=(["lat"],np.arange(-4,4.5,0.5)), lon=(["lon"],np.arange(1,7.5,0.5)))) # one year, 17 lat x 13 lon domain
    Filtered, Residual = meaneddy(Velocity, days=10, DataArray=True,timedim=["time"])

    
    INPUT:
       -> prop: 1, 2 or 3D array to filter
       -> days: number of days to set up the filter
       -> ndim: number of dimensions of the data [only used for DataArray=False, max:3]
       -> DataArray: True if prop is in xr.DataArray format
       -> dim: name of time dimension to filter (only used for DataArray=True)
   
    OUTPUT:
       -> m_prop: mean (filtered) part of the property
       -> p_prop: prime part of the property, essentially prop - m_prop
    
    v1 (February 2018)
    Cesar B. Rocha
    Dante C. Napolitano (dante.napolitano@legos.obs-mip.fr)
    
    v2 (December 2020)
    Dante C. Napolitano (dante.napolitano@legos.obs-mip.fr)
    Mariana M. Lage (mariana.lage@hereon.de)
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
