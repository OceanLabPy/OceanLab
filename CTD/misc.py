# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import seawater as sw
import matplotlib.pyplot as plt


def near(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return dat[idx][:how_many]

def argnear(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return idx
    
def argdistnear(x,y,xi,yi):
    idxs = []
    for xx,yy in zip(x,y):
        dists = np.sqrt((xi-xx)**2 + (yi-yy)**2)
        idxs.append(np.argmin(dists))
    return np.array(idxs)

def isopic_depth(DENS,PRES,isopic):
    '''
    This function looks for isopicnal depth from
    density vertical field based on pressure field.
    
    The method is based on linear interpolation only
    on vertical dimension.
    
    (2D-np.array,2D-np.array,float) -> 1D-np.array
    '''
    #checks if isopicnal value is on the range of data
    if (isopic>np.nanmax(DENS))|(isopic<np.nanmin(DENS)):
        #raise an error message
        raise ValueError('Isopicnal is out of range!')
    else:
        #defines a empty list
        p_ref = []
        #read each column
        for col in np.arange(DENS.shape[1]):
            #creates a interpolation function
            #this is based on density by pressure
            f = scint.interp1d(DENS[:,col],PRES[:,col])
            #returns the pressure based on isopicnal
            p_ref.append(f(isopic))
        #converts to numpy array
        p_ref = np.array(p_ref)
        return p_ref

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
    if type(df) <> pd.core.frame.DataFrame:
        #raise an Error
        raise ValueError('Data is not in the correct format (pandas.DataFrame)')              
    
    extrap_df = df.copy()
    
    nans_cnt = [0]
    #while the df has some nan value
    while True:
        nans_cnt.append(np.argwhere(pd.isnull(extrap_df)).shape[0])
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
            nans_cnt.append(np.argwhere(pd.isnull(extrap_df)).shape[0])
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
            new_values = ref_val.reshape(dif.shape)-dif
            #replace the nans by the new values
            extrap_df.iloc[np.hstack(nans),col] = new_values.reshape(ref_val.shape)
    return extrap_df




#
#sec = pd.read_pickle(u'/home/iury/TRABALHO/MARSEAL02/CTD/sections/SEAL02rad1')
#lat = sec.minor_xs('lat').mean().values
#lon = sec.minor_xs('lon').mean().values
#
#dist = np.hstack((0,np.cumsum(sw.dist(lat,lon)[0])))
#prof = sec.major_axis.values
#
#X,Z = np.meshgrid(dist,prof)
#
#gpan = extrap_all(sec.minor_xs('gpan'),lat=lat,lon=lon,inverse=True)
#
#gvel = sw.gvel(gpan,lat,lon)
#gvel_r = gvel-gvel[1100,:]
#
#X = X[:,:-1] + np.diff(X,axis=1)/2
#Z = Z[:,:-1]
#
#plt.ion()
#
#fig,(ax1,ax2) = plt.subplots(2,1)
#ax1.contourf(X,Z,gvel,np.arange(-2,2.1,0.1),cmap='seismic')
#ax1.set_ylim([0,1000])
#ax1.invert_yaxis()
#
#ax2.contourf(X,Z,gvel_r,levels=np.arange(-1,1.1,0.1),
#                cmap='seismic')
#ax2.set_ylim([0,1000])
#ax2.invert_yaxis()



#
#import pandas as pd
#import numpy as np
#df1 = pd.DataFrame(np.array([[1,1,1],[np.nan,3,2],[2,3,np.nan]]))
#df2 = extrap_all(df1.copy())
#
#X,Y =  np.meshgrid(df1.columns.values*1.,df1.index.values*1.)
#
#plt.figure()
#plt.contourf(X,Y,df1)
#plt.title('Not Extrapolated Data')
#plt.figure()
#plt.contourf(X,Y,df2)
#plt.title('Extrapolated Data')	
