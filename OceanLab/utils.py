import scipy.interpolate as scint
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import seawater as sw
import os
import pickle
from datetime import datetime
import requests
import shutil


# lonmin=-55.,lonmax=-22.,latmin=-36.,latmax=-4.6


def ETOPOget(lonmin,lonmax,latmin,latmax):
    url = "https://maps.ngdc.noaa.gov/mapviewer-support/wcs-proxy/wcs.groovy"+\
          "?filename=etopo1_bedrock.nc&request=getcoverage&version=1.0.0&"+\
          "service=wcs&coverage=etopo1_bedrock&CRS=EPSG:4326&format=netcdf&"+\
          "resx=0.016666666666666667&resy=0.016666666666666667&"+\
          "bbox={lonmin},{latmin},{lonmax},{latmax}"


    info = dict(lonmin=lonmin,lonmax=lonmax,latmin=latmin,latmax=latmax)

    r = requests.get(url.format(**info), stream=True)
    if r.status_code == 200:
        with open('tmp.nc', 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
        data = Dataset('tmp.nc')
        os.remove('tmp.nc')

    else:
        data = 'Internal Error'

    return data





#SAVING AND LOADING DATA IN PYTHON
def save_pickle(obj, name):
    '''
    Save python object as pickle binary.
    '''
    if name[-4:]=='.pkl':
        path = name
    else:
        path = name+'.pkl'

    with open(path, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pickle(name,enc='latin1'):
    '''
    Load python object as pickle binary.
    '''
    if name[-4:]=='.pkl':
        path = name
    else:
        path = name+'.pkl'
    with open(path, 'rb') as f:
        return pickle.load(f,encoding=enc)


def interp2_yaxis(X,Y,Z,yi,kind='linear'):
    '''
    This function interpolates a given data by second axis (Y-axis).

    Usage:
        interp2_axis1(2D-array,2D-array,2D-array,1D-array,string)
                                --> 2D-array,2D-array,2D-array

    '''
    #FUNCTION TO VERIFY IF (A/DE)SCENDENT OF Y AND YI IS CONSISTENT
    dif  = lambda dt: np.diff(dt[[0,-1]])[0]
    order = lambda dt: (dif(dt)/np.abs(dif(dt))).astype('int')
    #INTERPOLATION LAMBDA FUNCTION
    kw = {'kind':kind,'bounds_error':False}
    interp = lambda col: scint.interp1d(Y[:,col],Z[:,col],**kw)(yi)

    #INTERPOLATION
    Zi = np.vstack([interp(col) for col in np.arange(Z.shape[1])]).T
    #NEW GRID DATA
    Yi = np.array([yi,]*Zi.shape[1]).T
    Xi = np.array([X[0,:],]*yi.size)

    #CHANGING THE ORDER
    if order(yi)!=order(Y[:,0]):
        Zi = np.flipud(Zi); Yi = np.flipud(Yi)

    return Xi,Yi,Zi


def rsme(V_calc,V_obs):
    '''
    This function calculates the mean squared error.
    '''
    drange = lambda d: np.nanmax(d.ravel())-np.nanmin(d.ravel())

    n      = V_calc.size
    Vrange = drange(V_obs)

    err    = (np.sqrt(np.sum(np.subtract(V_calc,V_obs)**2)/n)/Vrange)*100
    return err


def nans(size):
    '''
    Create an array full of NaNs with some given shape.
    '''
    return np.zeros(size)*np.nan


def argdistnear(x,y,xi,yi):
    '''
    This function finds the index to nearest points in (xi,yi) from (x,y).

    '''
    idxs = [np.argmin(np.sqrt((xi-xx)**2 + (yi-yy)**2)) for xx,yy in zip(x,y)]
    idxs = np.array(idxs)
    return idxs


def download_bathy(lnd=-49,lnu=-33,ltd=-34,ltu=-20):
    '''
    This function downloads ETOPO1 data and make subset.

    lnd: lowest  longitude limit
    ltd:   "     latitude    "

    lnu: highest longitude limit
    ltu:    "    latitude    "
    '''
    #DECLARE URL PATH TO DOWNLOAD DATA
    arq = 'ETOPO1_Bed_g_gmt4.nc'
    url = 'http://www.ngdc.noaa.gov/thredds/dodsC/relief/ETOPO1/thredds/'+arq

    #READIND DATA
    etopodata = Dataset(url)

    #READING COORDINATES
    blon = (etopodata.variables['lon'][:])
    blat = etopodata.variables['lat'][:]

    #MAKING SUBSET
    condy = (blat>ltd)&(blat<ltu)
    condx = (blon>lnd)&(blon<lnu)
    topoin = etopodata.variables['z'][condy,:][:,condx]

    #MAKING GRID
    bLON,bLAT = np.meshgrid(blon[condx],blat[condy])

    return bLON,bLAT,topoin


#BOUNDARY CONDITION

def bathy_lims(PROF,lnu=-33,lnd=-49,
                ltu=-20,ltd=-34,step=1):

    bLON,bLAT,topoin = download_bathy(lnu=lnu,lnd=lnd,ltu=ltu,ltd=ltd)

    #PLOTTING
    plt.ioff()
    plt.figure()
    C = plt.contour(bLON,bLAT,-topoin,levels=[-PROF],colors='k')
    plt.close()

    #EXTRACTING COORDINATES
    C = np.vstack(C.allsegs[0])
    pseudo_lat = C[:,1]
    pseudo_lon = C[:,0]
    pseudo_dat = np.zeros(pseudo_lat.shape)

    return pseudo_lon[::step],pseudo_lat[::step],pseudo_dat[::step]


def deflagg(cast,flag=-9.990000e-29):
	'''
	pandas dataframe  -->  pandas dataframe
	Remove badflags data from pandas dataframe.

	flag default is -9.99e-29

	Example:
	deflagg(cast,flag=999)
	'''

	castdeflag = cast[cast.flag != flag]
	return castdeflag

def basename(fname):
    """Return filename without path.

    Examples
    --------
    >>> fname = '../test/data/FSI.txt.zip'
    >>> basename(fname)
    ('../test/data', 'FSI.txt', '.zip')
    """
    path, name = os.path.split(fname)
    name, ext = os.path.splitext(name)
    return path, name, ext

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)



def near(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return dat[idx][:how_many]

def argnear(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return idx[:how_many]


def select_rad(pts,lon,lat):
    args=np.array([])
    for est in pts:
        c=[]
        arg=np.array([0])
        for i in np.arange(0,lon.size):
            c.append(sw.dist([lat[i],est[1]],[lon[i],est[0]]))
        arg=np.argmin(np.array(c).T[0][0])
        args=np.append(args,arg)
    return args

###############################################################################
# EXTRAPOLATION
###############################################################################

def extrap_all(df,lat=[],lon=[],inverse=True,wgt=0.5):
    '''
    This function extrapolate a section of some property
    organized as a pandas.DataFrame using the horizontal
    gradient of attribute. 'inverse' is the boolean variable
    that flag if need to extrapolate in 2 ways of horizontal
    data. 'wgt' is the weigth applied to gradient.

    (pandas.DataFrame) => pandas.DataFrame

    Give a try:

    import pandas as pd
    import numpy as np
    import seawater as sw
    import matplotlib.pyplot as plt

    rand = np.random.randint(15,30,30)*1.
    rand[np.random.randint(0,rand.size,13)] = np.nan
    rand = rand.reshape((6,5))
    data = pd.DataFrame(rand,columns=list('ABCDE'))

    ficlat,ficlon = np.array([-10,-10,-10,-10,-10]),np.arange(-38,-33)

    data_ext = extrap_all(data,lat=ficlat,lon=ficlon)

    dist = np.hstack([0,np.cumsum(sw.dist(ficlat,ficlon)[0])])

    D,P = np.meshgrid(dist,-data.index.values)

    fig,(a1,a2) = plt.subplots(2,1)
    a1.contourf(D,P,data)
    a1.set_title('Before')
    a2.contourf(D,P,data_ext)
    a2.set_title('After')
    plt.show()


    '''
    #if the df is not a pandas.DataFrame
    if type(df) != pd.core.frame.DataFrame:
        #raise an Error
        raise ValueError('Data is not in the correct format (pandas.DataFrame)')

    extrap_df = df.copy()

    nans_cnt = [0]
    #while the df has some nan value
    while True:
        nans_cnt.append(np.argwhere(pd.isnull(extrap_df).values).shape[0])
        if nans_cnt[-1]!=nans_cnt[-2]:
            extrap_df = extrap_gradient(extrap_df,lat,lon,wgt=wgt)
        else:
            break

    if inverse:
        nans_cnt = [0]
        extrap_df = extrap_df.iloc[:,::-1]
        lat = lat[::-1]
        lon = lon[::-1]
        #while the df has some nan value
        while True:
            nans_cnt.append(np.argwhere(pd.isnull(extrap_df).values).shape[0])
            if nans_cnt[-1]!=nans_cnt[-2]:
                extrap_df = extrap_gradient(extrap_df,lat,lon,wgt=wgt)
            else:
                break

        extrap_df = extrap_df.iloc[:,::-1]
        lat = lat[::-1]
        lon = lon[::-1]

    return extrap_df


def extrap_gradient(df,lat=[],lon=[],wgt=0.5):

    ##if the df is not a pandas.DataFrame
    #if type(df) <> pd.core.frame.DataFrame:
    #    #raise an Error
    #    raise ValueError('Data is not in the correct format (pandas.DataFrame)')

    #calculate the distance vector from the first
    dist = np.hstack((0,np.cumsum(sw.dist(lat,lon)[0])))
    #calculate the distance between the points
    dxs = sw.dist(lat,lon)[0]


    #calculate the horizontal gradient between the points
    gradient = ((np.diff(df))/dxs)*wgt
    # *1. is to make sure we
    #are working with float numbers

    extrap_df = df.copy()
    #Assuming the points with more nan values is the first
    #we just read until the second from last column
    for col in np.arange(extrap_df.shape[1]-2):
        #find where is the nan values along the column
        nans = np.argwhere(np.isnan(extrap_df.iloc[:,col]))

        if nans.tolist()!=[]:
            #the gradient per unit of distance plus the distance
            #between the column and the next one
            dif = gradient[nans,col+1]*dxs[col]
            #find the next value to apply the coef
            ref_val = extrap_df.iloc[np.hstack(nans),col+1]
            #calculate the new value to the nans in col based on values in col+1
            #new_values = ref_val.reshape(dif.shape) - dif
            new_values = ref_val - dif.squeeze()
            #replace the nans by the new values
            #extrap_df.iloc[np.hstack(nans),col] = new_values.reshape(ref_val.shape)
            extrap_df.iloc[np.hstack(nans),col] = new_values
    return extrap_df
