# -*- coding: utf-8 -*-
import os
import pandas as pd
from glob import glob
import numpy as np
import gsw
import seawater as sw
from collections import OrderedDict

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

