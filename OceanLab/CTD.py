# -*- coding: utf-8 -*-
import os
import pandas as pd
from glob import glob
import numpy as np
import gsw
import seawater as sw
from collections import OrderedDict
import scipy.interpolate as scint

#DYNAMICAL MODES

def eqmodes(N2,z,nm,pmodes=False):
    '''
    This function computes the equatorial velocity modes
    
    ========================================================
    
    Input:
        
        N2 - Brunt-Vaisala frequency data array  
    
        z - Depth data array (equaly spaced)
    
        lat - Latitude scalar 
    
        nm - Number of modes to be computed
    
        pmodes - If the return of pressure modes is required. 
                 Default is False
    
    ========================================================

    Output:
        
        Si - Equatorial modes matrix with MxN dimension being:
                M = z array size
                N = nm

        Rdi - Deformation Radii array
        
        Fi - Pressure modes matrix with MxN dimension being:
                M = z array size
                N = nm
                
            Returned only if input pmodes=True
    
    made by Hélio Almeida, Iury Sousa and Wandrey Watanabe
    Laboratório de Dinâmica Oceânica - Universidade de São Paulo
                                2016
    '''

    #needed to the problem
    lat = 0
    #nm will be the number of baroclinic modes
    nm -= 1 #Barotropic mode will be added

    #defines the orthonormalization function
    onorm      = lambda f: f/np.sqrt(dz*((np.array(f[1:])**2+\
                                    np.array(f[:-1])**2)/(2*H)).sum())
    # defines function to keep consistency multiply by plus/minus 
    # to make the leading term positive
    plus_minus = lambda f: -f if (np.sign(f[1]) < 0) else f
    		
    		
    dz = np.abs(z[1] - z[0])

    # assembling matrices
    # N2 matrix
    N2 = np.diag(N2[1:-1],0)
    # 2nd order difference matrix    
    A  = np.diag(np.ones(z.size-2)*-2.,0)+\
         np.diag(np.ones(z.size-3)*1.,-1)+\
         np.diag(np.ones(z.size-3)*1.,1)  
    A  = A/(dz**2)

    A  = np.matrix(A)
    N2 = np.matrix(N2)
    
    #C  = A*N2
    N02 = -1/N2
    N02[np.isinf(N02)]=0
    C  = N02*A
    
    # solve the eigenvalue problem
    egval,egvec = np.linalg.eig(C)

    i     = np.argsort(egval)
    egval = egval[i]
    egvec = egvec[:,i]
    
    # truncate the eigenvalues and eigenvectors for the number of modes needed
    ei = egval[:nm]
    Si = egvec[:,:nm]
    
    # Applying Dirichlet boundary condition at top and bottom
    # adding a row of zeros at bottom and top
    Si   = np.append(np.matrix(np.zeros(nm)),Si, axis=0)
    Si   = np.append(Si,np.matrix(np.zeros(nm)), axis=0)
    Si   = np.array(Si)


    H  = np.abs(z).max()
    # normalizing to get orthonormal modes 
    Si = np.array(map(onorm,Si.T)).T
    # to keep consistency multiply by plus/minus 
    # to make the leading term positive
    Si = np.array(map(plus_minus,Si.T)).T      

    # compute the deformation radii [km]    
    beta = (7.2921150e-5*2*np.cos(np.deg2rad(lat*1.)))/6371000
    c    = np.sqrt(1/ei)
    radii= np.sqrt(c/beta)
    #append external deformation radius
    radii= np.hstack([(np.sqrt(9.81*H)/np.abs(sw.f(lat))),radii])
    #converting to km
    radii*=1e-3
    
    
    #BAROTROPIC MODE
    no = 1
    fb = np.ones((Si.shape[0],1))*no
    sb = np.expand_dims(np.linspace(0,no+1,Si.shape[0]), axis=1)
    
    # trying to compute the pmodes based on the polarization
    # relation between velocity/pressure modes
    #
    # F_j(z)=-g.he_j d/dzS_j
    #
    ####
    if pmodes==True:
        Fi=np.zeros(Si.shape)
        for i in np.arange(nm):
            Fi[1:-1,i] =\
            (-1/ei[i])*((Si[1:-1,i]-Si[0:-2,i])/dz)   

        #Aplying Neuman boundary condition d/dzFi=o @0,-H
        Fi[0],Fi[-1]=Fi[1],Fi[-2]
        

        # normalizing to get orthonormal modes 
        Fi = np.array(map(onorm,Fi.T)).T
        # to keep consistency multiply by plus/minus 
        # to make the leading term positive
        Fi = np.array(map(plus_minus,Fi.T)).T
        
        Si = np.hstack([sb,Si])
        Fi = np.hstack([fb,Fi])

        return Si,radii,Fi 
    else:
        Si = np.hstack([sb,Si])  
        return Si,radii   


#MISC

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    
    
#EXTRACTION

def get_lonlatnames(fname,lathint='Latitude =',lonhint='Longitude =',
                          lonline=[],latline=[]):
    f = open(fname)
    header, config, names = [], [], []
    for k, line in enumerate(f.readlines()):
        line = line.strip()
        if line.startswith('*'):  # Get header.
            header.append(line)
        if line.startswith('#'):  # Get configuration file.
            config.append(line)
        if '# name' in line:
            names.append(line.split('=')[-1].split(':')[0][1:])
        if (lonline==[])|(latline==[]):
            latcond = lathint in line
            loncond = lonhint in line
        else:
            latcond = k==latline
            loncond = k==lonline

        if latcond:
            hemisphere = line.split()[-1]
            deg = float(line.split()[-3])
            minute = float(line.split()[-2])
            lat = (deg + minute / 60.)
            if hemisphere == 'S':
                lat *= -1
            else:
                raise ValueError("Latitude not recognized.")

        if loncond:
            hemisphere = line.split()[-1]
            deg = float(line.split()[-3])
            minute = float(line.split()[-2])
            lon = (deg + minute / 60.)
            if hemisphere == 'W':
                lon *= -1
            else:
                raise ValueError("Latitude not recognized.")
        if line == '*END*':  # Get end of header.
            skiprows = k + 1
            break
    return lon,lat,names,skiprows



#PROCESSING  
    

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


def loopedit(cast):
	'''
	pandas dataframe  -->  pandas dataframe
	Remove loop data from pandas CTD dataframe.
	This loopedit function delete only the data
	from returnind CTD.

	Methods:
	    'press' = remove data fro negative np.diff(press)
	    'dz/dtM' = remove data from descent rate 

	
	'''


	try:
	    flag = cast['dz/dtM'].values>0
	    castedited = cast.iloc[flag,:]
	except:
	    flag = np.hstack([1,np.diff(cast.index.values)])>0
	    castedited = cast.iloc[flag,:]
	    
	return castedited
	
	


def abv_water(cast,maxprof=11000):
	'''
	returns down and up cast
	'''
	castdig = cast[cast.index.values<11000]
	castdig = castdig[castdig.index.values>0]
	down = cast.iloc[:cast.index.argmax()]
	up = cast.iloc[cast.index.argmax():][::-1]
	return down,up




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

def ctdread(fname,press_name='prDM',lathint='Latitude =',
                lonhint='Longitude =',down_cast=True,
                latline=[],lonline=[]):

	
	lon,lat,names,skiprows = get_lonlatnames(fname,lathint=lathint,
	                       lonhint=lonhint,lonline=lonline,latline=latline)


	cast = pd.read_csv(fname,skiprows=skiprows,
			names=names,delim_whitespace=True)

	cast.set_index(press_name,drop=True,inplace=True)
	cast.index.name = "Pressure [db]"
	cast = deflagg(cast)
	dwn,up = abv_water(cast)
	
	if down_cast:
	   return lon,lat,dwn
	else:
	   return lon,lat,up

def despike(self,propname,block,wgth=2):
	prop = np.array(self[propname])
	wint = rolling_window(prop,block)
	stdt = wgth * wint.std(axis=1)
	meant = wint.mean(axis=1)
	stdt = np.hstack([np.tile(stdt[0], (block - 1)/2),stdt,np.tile(stdt[-1], (block - 1)/2)])
	meant = np.hstack([np.tile(meant[0], (block - 1)/2),meant,np.tile(meant[-1], (block - 1)/2)])
	self = self[np.abs(self[propname]-meant)<stdt]
	return self

def hann_filter(self,propname,block):
    '''
    This function apply Hanning Window filter to some item
    named 'propname' 
    
    '''
    def hann(x):
        return (x*np.hanning(x.size)).sum()/np.hanning(x.size).sum()
        
    filtered_na = pd.rolling_apply(self[propname],block,hann,center=True)
    #Fill the head and tail of values that does not got the filter
    self[propname] = filtered_na.fillna(self[propname])
    
    return self
    

def binning(self,delta=1.):
        start = np.floor(self.index[0])
        end = np.ceil(self.index[-1])
        shift = delta / 2.  # To get centered bins.
	bins = np.arange(start,end,1.)-shift
	binned = self.groupby(np.digitize(self.index.values.astype('float'),bins)).mean()
	return binned

def ctdproc(lista,temp_name='t068C',
        lathint='Latitude =',lonhint='Longitude =',
        cond_name='c0S/m',press_name='prDM',down_cast=True,
        looped=True,hann_f=False,hann_block=20,hann_times=2,
        latline=[],lonline=[]):
    '''
    This function do the basic proccess to all .cnv CTD data from
    given list.
    '''
    for fname in lista:
        
   	lon,lat,data = ctdread(fname,press_name=press_name,
   	                down_cast=down_cast,lathint=lathint,
   	                lonhint=lonhint,lonline=lonline,latline=latline)
   	
   	if looped:
           	data = loopedit(data)
           	
   	dataname = basename(fname)[1]
    
   	if (data.shape[0]<101)&(data.shape[0]>10): # se o tamanho do perfil for com menos de 101 medidas
    
  		if (data.shape[0]/2)%2 == 0: # caso a metade dos dados seja par
 			blk = (data.shape[0]/2)+1 # bloco = a metade +1
  		else:
 			blk = data.shape[0]/2 # se for impar o bloco e a metade
    
  		# remove spikes dos perfis de temperatura e condutividade
  		data = despike(data,propname=temp_name,block=blk,wgth=2)
  		data = despike(data,propname=cond_name,block=blk,wgth=2)
   	elif data.shape[0]>=101:
  		# para perfis com mais de 101 medidas, utiliza-se blocos de 101
  		data = despike(data,propname=temp_name,block=101,wgth=2)
  		data = despike(data,propname=cond_name,block=101,wgth=2)
  	else:
  	         print 'radial muito rasa'

   	# realiza média em caixa de 1 metro
   	data = binning(data,delta=1.)
   	if temp_name=='t068C':
   	    data['t090C'] = gsw.t90_from_t68(data['t068C'])
   	       	    
   	data['sp'] = gsw.SP_from_C(data[cond_name]*10,data['t090C'],data.index.values)

   	if hann_f:
   	    times=0
   	    while times<hann_times:
           	    data = hann_filter(data,'t090C',hann_block)
           	    data = hann_filter(data,'sp',hann_block)
           	    times +=1

   	data['pt'] = sw.ptmp(data['sp'],data['t090C'],data.index.values)
   	#data['ct'] = gsw.CT_from_pt(data['sa'],data['pt'])
   	data['psigma0'] = sw.pden(data['sp'],data['t090C'],data.index.values,pr=0)-1000
   	data['psigma1'] = sw.pden(data['sp'],data['t090C'],data.index.values,pr=1000)-1000
   	data['psigma2'] = sw.pden(data['sp'],data['t090C'],data.index.values,pr=2000)-1000
   	data['gpan'] = sw.gpan(data['sp'],data['t090C'],data.index.values)
   	data['lat'] = lat
   	data['lon'] = lon
   	
   	data.to_pickle(os.path.split(fname)[0]+'/'+os.path.splitext(os.path.split(fname)[1])[0])
   	
   	print dataname


def return_section(directory,ext='*.cnv'):
	'''
	This function reads all the files from directory
	that has no extension as a pandas pickle, sorted
	by name and return a panel section and its lat 
	and lon.
	'''

	# this is the list of cnv files from directory
	# sorted by name
	stas = np.sort(glob(os.path.join(directory,ext)))
	# define section as an empty ordered dictionary
	section = OrderedDict()
	# starts the loop
	for pf in stas:
           	# read the file with the same name and no extension
           	# because pandas pickle has no extensio
		data = pd.read_pickle(pf.split('.')[0])
		
		#if name==None:
  #                 	# define the name of the file for the dict
  #                 	# again, the name has no extension (use split)
          	name = os.path.basename(pf).split('.')[0]
           	# give to section an update, so this dict
           	# will have a variable with the same file name
		section.update({name: data})
		
	# transform from dict to panel
	section = pd.Panel.fromDict(section)
	
	section = section.fillna(method='backfill') 
	
	# return the lat, lon and section
	return section


### CALC SECTION ##############################################################

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


def isopic_depth(DENS,PRES,isopic,med=False):
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
        if med:
            #calculates the medium point
            p_ref = p_ref[:-1]+np.diff(p_ref)/2
        #converts to numpy array as rounded index
        p_ref = np.array(p_ref).round().astype('int')
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
###############################################################################

#		
#lista = glob('/home/iury/TRABALHO/MARSEAL01/CTD_hann/*.cnv')
#
#ctdproc(lista,temp_name='t090C',looped=True,hann_block=7,hann_times=2,hann_f=True)
#
#lista = glob('/home/iury/TRABALHO/MARSEAL01/CTD/*.cnv')
#
#ctdproc(lista,temp_name='t090C',looped=True,hann_f=False)
#
#
#lista = glob('/home/iury/TRABALHO/MARSEAL02/CTD_hann/*.cnv')
#
#ctdproc(lista,temp_name='t090C',looped=True,hann_block=7,hann_times=2,hann_f=True)
#
#lista = glob('/home/iury/TRABALHO/MARSEAL02/CTD/*.cnv')
#
#ctdproc(lista,temp_name='t090C',looped=True,hann_f=False)


#
#
#sec = return_section('/home/iury/TRABALHO/MARSEAL01/CTD/')
#sec2 = return_section('/home/iury/TRABALHO/MARSEAL01/CTD_hann/')
#
#plt.plot(sec.minor_xs('pt'),'b')
#plt.plot(sec2.minor_xs('pt'),'r')


#lon,lat,data = ctdread(fname,press_name='prDM',latline=14,lonline=15)

#path = '/home/iury/Copy/TCC/Rotinas/dados'

#lista = [os.path.join(dirpath, f)
#    for dirpath, dirnames, files in os.walk(path)
#    for f in fnmatch.filter(files, '*.cnv')]
#
#
#path = '/home/iury/TRABALHO/POTIGUAR/CTD/'
#
#lista = np.sort(glob(path+'*.cnv'))
#
#ctdproc(lista,latline=14,lonline=15,loopedit=False)
#
#
#sec = return_section('/home/iury/TRABALHO/POTIGUAR/CTD/')

