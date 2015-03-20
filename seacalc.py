# -*- coding: utf-8 -*-
#Most of these functions is used to calculate geostrophic velocity, interpolation,
#extrapolation methods, calculate stream function, sverdrup transport, etc...
#Made by Iury Sousa - São Paulo/Brazil

#Some of these functions were based or inspired by other authors
#like Filipe Fernandes and his scripts published on GitHub or his 
#website https://ocefpaf.github.io/python4oceanographers/
#I have a lot to thank for his help, manly giving his time to
#show his scripts on web.

import os
from glob import glob

import pandas as pd
from pandas import read_table,Series,Index
#from pandas import Panel

import numpy as np
import numpy.ma as ma
from collections import OrderedDict

#import scipy.interpolate as interp
import scipy.io as sio
import scipy.interpolate as scint

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
from matplotlib import gridspec
from matplotlib.transforms import Bbox
from matplotlib.colors import BoundaryNorm
#import fnmatch

from math import atan2, degrees, pi


import gsw
import seawater as sw



def vectoa(xc,yc,x,y,u,v,corrlen,err,b):
	#%function psi=vectoa(xc,yc,x,y,u,v,corrlen,err)
	#%       (xc,yc) are vectors (row or column) of interpolation points
	#%       (x,y)   are vectors (row or column) of observation points
	#%       (u,v)   are matrices of east and north components with each day
	#%               entered column-wise or row-wise
	#%       corrlen,err are correlation length scales and error for a
	#%               gaussian streamfunction covariance function
	#%       psi     is streamfunction at    (xc,yc) returned in same
	#%               format as u or v
	#PYTHON VERSION MADE BY Iury Simões - Oceanography from Universidade Federal do Ceará
	# 2014

	n = len(x)
	#mu,nu=u.shape
	#if (mu==n):
	#        u=np.vstack((u,v))        #data entered column-wise
	#elif(nu==n):
	u=np.array([np.hstack((u,v))]).T    #data entered row-wise
	#else:
	#        print('Data array has incorrect size, does not match (x,y).')
	#	pass

	#% angles and distances
	pp = -np.tile(y,(n,1)).T+np.tile(y,(n,1)),-np.tile(x,(n,1)).T+np.tile(x,(n,1)) 
	t = []
	for ii,jj in zip(pp[0].ravel(),pp[1].ravel()):
		t.append(atan2(ii,jj))
	t = np.array(t)
	t.shape = pp[0].shape	
	d2=((np.tile(x,(n,1)).T-np.tile(x,(n,1)))**2+(np.tile(y,(n,1)).T-np.tile(y,(n,1)))**2)

        

	lambd = 1/(corrlen**2)
	bmo=b*err/lambd
	R=np.exp(-lambd*d2)        #%longitudinal
	S=R*(1-2*lambd*d2)+bmo #%transverse
	R=R+bmo

	A=np.zeros((2*n,2*n))

	A[0:n,0:n]=(np.cos(t)**2)*(R-S)+S
	A[0:n,n:2*n]=np.cos(t)*np.sin(t)*(R-S)
	A[n:2*n,0:n]=A[0:n,n:2*n]
	A[n:2*n,n:2*n]=(np.sin(t)**2)*(R-S)+S
	A=A+err*np.eye(A.shape[0])

	#% angles and distances
	nv1,nv2 =xc.shape
	nv=nv1*nv2

	
	xc = xc.T.ravel()
	yc = yc.T.ravel()

	
	#% angles and distances
	ppc = -np.tile(yc,(n,1)).T+np.tile(y,(nv,1)),-np.tile(xc,(n,1)).T+np.tile(x,(nv,1)) 
	tc = []
	for ii,jj in zip(ppc[0].ravel(),ppc[1].ravel()):
		tc.append(atan2(ii,jj))
	tc = np.array(tc)
	tc.shape = ppc[0].shape	
	d2=((np.tile(xc,(n,1)).T-np.tile(x,(nv,1)))**2+(np.tile(yc,(n,1)).T-np.tile(y,(nv,1)))**2)
	R=np.exp(-lambd*d2)+bmo;
	
	
	P=np.zeros((nv,2*n))
	#%streamfunction-velocity covariance
	P[:,0:n]=np.sin(tc)*np.sqrt(d2)*R;
	P[:,n:2*n]=-np.cos(tc)*np.sqrt(d2)*R;

	psi=np.dot(P,np.linalg.solve(A,u))   #%uses this line if u is full
	#%adjust if column-wise
#	if(nu==n):       
#	        psi=psi.T
	return psi        
        


	


def closest_node(node, nodes):
    '''
    This function find the index of nearest point considering
    2 component coordinates.
    
    (2-element list or array,more-than-2-elements list or array) => number
    '''
    #convert to array
    nodes = np.asarray(nodes)
    #calculate the sum of quadratic difference between
    #points of each coordinate
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)
    




def extrap_gradient(df,coords=False,lat=[],lon=[]):
	'''
	This function extrapolate a section of some property
	organized as a pandas.DataFrame using the horizontal
	gradient of attribute. coords is the boolean variable
	that flag if the DataFrame columns is the distance vector
	between the points or if the latitudes and longitudes
	is given.
	
	(pandas.DataFrame) => pandas.DataFrame
	
	Give a try:
	import pandas as pd
	import numpy as np
	df1 = pd.DataFrame(np.array([[1,1,1],[np.nan,3,2],[2,3,np.nan]]),columns=np.array([0,100,600]))
	df2 = extrap_gradient(df1.copy())
	
	X,Y =  np.meshgrid(df1.columns.values*1.,df1.index.values*1.)
	
	plt.figure()
	plt.contourf(X,Y,df1)
	plt.title('Not Extrapolated Data')
	plt.figure()
	plt.contourf(X,Y,df2)
	plt.title('Extrapolated Data')	

	
	'''
	#if the df is not a pandas.DataFrame
	if type(df) <> pd.core.frame.DataFrame:
		#raise an Error
		raise ValueError('Data is not in the correct format (pandas.DataFrame)')
	
	#if the coordinates is given as latitudes and longitudes
	if coords:
		#calculate the distance vector from the first
		dist = np.hstack((0,np.cumsum(sw.dist(lat,lon)[0])))
		#calculate the distance between the points
		dxs = sw.dist(lat,lon)[0]
	else:
		#other cases
		#define the distance as the columns of DataFrame
		dist = df.columns.values
		#and the distance between the points as its diff
		dxs = np.diff(dist)

	#begin a counter
	cnt = 1
	#while the df has some nan value
	while np.any(np.isnan(df)):
		#calculate the horizontal gradient between the points
		gradient = ((np.diff(df.values))/dxs)*1. # *1. is to make sure we 
									      #are working with float numbers
	
	#Assuming the points with more nan values is the first
	#we just read until the second from last column
		for col in np.arange(df.shape[1]-2):
			#find where is the nan values along the column
			nans = np.argwhere(np.isnan(df.iloc[:,col]))
			#inf there's no nan value
			if nans.size == 0:
				#print a warning
				print 'Not nans in this column'
			else:
				#the gradient per unit of distance plus the distance
				#between the column and the next one
				dif = gradient[nans,col+1]*dxs[col]
				#find the next value to apply the coef
				ref_val = df.iloc[np.hstack(nans),col+1]
				#calculate the new value to the nans in col based on values in col+1
				new_values = ref_val.reshape(dif.shape)-dif
				#replace the nans by the new values
				df.iloc[np.hstack(nans),col] = new_values.reshape(ref_val.shape)
		#try one more time
		cnt +=1
		
		#Even tehre's NaN values after 10 times, the method change to
		#the nearest value
		if cnt > 10:
			print 'Radial com estação única com profundidade máxima'
			df = df.interpolate(method='nearest', axis=0).ffill().bfill()
			break #break the while loop
			
	return df
		


def interpgrid(self,dya=None,interp=False,d4=None):

	a = self
	
	xa = a.columns.values.astype('float')
	
	ya = a.index.values.astype('float')

	X,Y = np.meshgrid(xa,ya)
	
	if interp:
		xa2 = d4
	else:
		xa2 = xa
	
	try: 
		Xn,Yn = np.meshgrid(xa2,np.arange(0,dya+1))
	except:
		Xn,Yn = X,Y
		
	a2 = a.values

	mask = ~np.isnan(a2.ravel())
	points = zip(X.ravel()[mask],Y.ravel()[mask])
	new_points = zip(Xn.ravel(),Yn.ravel())

	a3 = scint.griddata(points,a2.ravel()[mask],new_points)
	
	if interp:
		a3.shape = (dya+1,xa2.size)
		a3 = pd.DataFrame(a3,index=np.arange(0,dya+1),columns=xa2)
			
	else:
		try:
			a3.shape = (dya+1,a2.shape[1])
			a3 = pd.DataFrame(a3,index=np.arange(0,dya+1),columns=xa2)		
		except:
			a3.shape = a2.shape
			a3 = pd.DataFrame(a3,index=ya,columns=xa2)
	

	return a3
	
	
	
def mdr_fixed(vsec_data,adcp_data,refprof):
	reff = adcp_data.iloc[refprof,:].values-vsec_data.iloc[refprof,:].values
	vsec_data += reff
	return vsec_data

def mdr(vsec,adcp):
#	vsec,lat,lon,section = extdat(rad,sup=True)
#	lat_adcp,lon_adcp,adcp = adcpvel(lat,lon,'/home/iury/Copy/TCC/new_contour_one',dya=600)

	adcprange = np.nanmax(np.nanmax(adcp))-np.nanmin(np.nanmin(adcp))


	pref = []
	for j in np.arange(0,vsec.shape[1]):
		p_ini = 100
		minrquem = []
		for i in np.arange(p_ini,250):
			vsecref = mdr_fixed(vsec,adcp,i)
			rquem = (np.sqrt(np.subtract(vsecref.iloc[0:601,j].values,adcp.iloc[:,j].values)**2)/adcprange)*100
			minrquem.append(np.nanmean(rquem))
		
		if np.all(np.isnan(minrquem)):
			pref.append(600-p_ini)
		else:
			pref.append(np.nanargmin(minrquem))
		
	pref = np.array(pref)		
	pref += p_ini

	velref = []
	for i,j in zip(pref.astype('int'),np.arange(0,adcp.shape[1])):
		reff = adcp.iloc[i,j]-vsec.iloc[i,j]
		velref.append(reff)
		
	vsec += velref
	return vsec,pref	
	

def doublesum(self):
	'''
	This function calculates a double sum of some data.
	is used to avoid write always numpy.sum(numpy.sum(ARRAY))
	'''
	summ = np.sum(np.sum(self))
	return summ

def sverdrup(section,secvel,radcomplicada=False):
	
	lat = section.minor_xs('lat').mean().values
	lon = section.minor_xs('lon').mean().values
	
	lonv = secvel.lon.mean().values
	latv = secvel.lat.mean().values

	if radcomplicada:
		latv[3] = -2.0933934592636545
		lonv[3] = -38.442823159524877

	dist = sw.dist(lat,lon)[0]
	dist = dist*1000

	a,bathy = etopo_sec(latv,lonv)
	bathy = -np.array(bathy).astype('int')
	
	
	dens0 = secvel.dens0
	dens1 = secvel.dens1

	
	a,pres = np.meshgrid(a,dens0.index.values)
	pres = pd.DataFrame(pres,columns=dens0.columns.values)
	batmask = pres<bathy
	

	mask1 = (dens0.values<24.5)&batmask
	mask2 = ((dens0.values>24.5)&(dens0.values<26.8))&batmask
	mask3 = ((dens0.values>26.8)&(dens1.values<32.15))&batmask
	mask4 = (dens1.values>32.15)&batmask
	

	svadcp = secvel.adcp*dist
	svmdr = secvel.velmdr*dist
#	sviso = secvel.veliso*dist
	
	tab = OrderedDict()
	
	pos,neg = [],[]
	
	pos = np.array([doublesum(svadcp[svadcp>0][mask1]),
	doublesum(svadcp[svadcp>0][mask2]),
	doublesum(svadcp[svadcp>0][mask3]),
	doublesum(svadcp[svadcp>0][mask4])])/10**6
	
	neg = np.array([doublesum(svadcp[svadcp<0][mask1]),
	doublesum(svadcp[svadcp<0][mask2]),
	doublesum(svadcp[svadcp<0][mask3]),
	doublesum(svadcp[svadcp<0][mask4])])/10**6
	
	tab.update({'adcp':pd.DataFrame(np.vstack((pos,neg)).T,
	index=['upper','ACAS','AIA','APAN'],
	columns=['pos','neg'])})
	
	pos,neg = [],[]
	
	pos = np.array([doublesum(svmdr[svmdr>0][mask1]),
	doublesum(svmdr[svmdr>0][mask2]),
	doublesum(svmdr[svmdr>0][mask3]),
	doublesum(svmdr[svmdr>0][mask4])])/10**6
	
	neg = np.array([doublesum(svmdr[svmdr<0][mask1]),
	doublesum(svmdr[svmdr<0][mask2]),
	doublesum(svmdr[svmdr<0][mask3]),
	doublesum(svmdr[svmdr<0][mask4])])/10**6
	
	tab.update({'mdr':pd.DataFrame(np.vstack((pos,neg)).T,
	index=['upper','ACAS','AIA','APAN'],
	columns=['pos','neg'])})
	
#	pos,neg = [],[]
#	
#	pos = np.array([doublesum(sviso[sviso>0][mask1]),
#	doublesum(sviso[sviso>0][mask2]),
#	doublesum(sviso[sviso>0][mask3]),
#	doublesum(sviso[sviso>0][mask4])])/10**6
#	
#	neg = np.array([doublesum(sviso[sviso<0][mask1]),
#	doublesum(sviso[sviso<0][mask2]),
#	doublesum(sviso[sviso<0][mask3]),
#	doublesum(sviso[sviso<0][mask4])])/10**6
#	
#	tab.update({'iso':pd.DataFrame(np.vstack((pos,neg)).T,
#	index=['upper','ACAS','AIA','APAN'],
#	columns=['pos','neg'])})

	tab = pd.Panel.from_dict(tab)
	return tab
	
	
def adcpvel(lat,lon,filepath,dya=None,interp=False):

	x,y,dday,u,v = [],[],[],[],[]
	xref,ang = [],[]
	#lendo dados de saída do adcp		
	uv = sio.loadmat(os.path.join(filepath,'contour_uv.mat'))
	xy = sio.loadmat(os.path.join(filepath,'contour_xy.mat'))
	
	#definindo longitute (-360) e latitude
	x= xy['xyt'][0,:]-360.
	y= xy['xyt'][1,:]
	#definindo variavel tempo
	dday= xy['xyt'][2,:]
	#definindo u e v
	u=uv['uv'][:,::2]
	v=uv['uv'][:,1::2]
	
	#retirando valores nan de lat e lon
	y = y[~np.isnan(x)]
	dday = dday[~np.isnan(x)]
	u = u[:,~np.isnan(x)]
	v = v[:,~np.isnan(x)]
	x = x[~np.isnan(x)]
	
	mask = np.all(np.isnan(u),axis=0)
	
	y = y[~mask]
	dday = dday[~mask]
	u = u[:,~mask]
	v = v[:,~mask]
	x = x[~mask]
	
	lats,lons = np.array([]),np.array([])
	U,V = np.array([]),np.array([])

	try: 
		leng = lat.size-1
	except:
		leng = np.array(lat).size-1
	

#	FINDING THE CLOSEST POINTS TO EXTREMES OF RADIAL	
	nodes = np.array(zip(x,y))
	t1 = closest_node([lon[0],lat[0]],nodes)
	t2 = closest_node([lon[-1],lat[-1]],nodes)

#	IDENTIFYING IF THE ADCP DATA IS INVERTED
	if t1>t2:
		t1,t2 = t2,t1
		invert=True
	else:
		invert=False


	latv = lat[:-1]+np.diff(lat)/2
	lonv = lon[:-1]+np.diff(lon)/2

#	CALCULATE THE DISTANCE AND AVERAGE ANGLE BETWEEN CTD STATIONS
	dist,ang = sw.dist(lat,lon)
	dist = np.hstack([0,np.cumsum(dist)])	
	ang = np.mean(np.radians(ang))


#	PRINT THE AVERAGE ANGLE
	print 'Ângulo = %1.2f' %(np.degrees(ang))


#	SELECTING INSIDE RADIAL DATA		
	yy = y[t1:t2]
	xx = x[t1:t2]
	
	uu = u[:,t1:t2]
	vv = v[:,t1:t2]

	if invert:
		yy = np.flipud(yy)
		xx = np.flipud(xx)
		uu = np.fliplr(uu)
		vv = np.fliplr(vv)	
	
#	CALCULATE THE DISTANCE OF EACH POINT TO THE ORIGIN
	dd,adcpang = sw.dist(np.array(zip(np.tile(lat[0],yy.shape),yy)).ravel(),
		     np.array(zip(np.tile(lon[0],xx.shape),xx)).ravel())
	dd = dd[::2]
	adcpang = np.radians(adcpang[::2])

	dd = dd*np.cos(adcpang-ang)
	
	
	uu = pd.DataFrame(uu,index=np.hstack(xy['zc']).astype('int'),columns=dd)
	vv = pd.DataFrame(vv,index=np.hstack(xy['zc']).astype('int'),columns=dd)
	
	plt.figure()
	plt.plot(x,y)
	plt.plot(xx,yy,linewidth=0,marker='o',markersize=3)
	
	
	if interp:
		
		distv = dist[:-1]+np.diff(dist)/2
		uu = interpgrid(uu,dya=dya,interp=True,d4=distv)
		vv = interpgrid(vv,dya=dya,interp=True,d4=distv)
		
		uu = interpgrid(uu,dya=dya)
		vv = interpgrid(vv,dya=dya)
					
	else:
		distv=dd
		uu = interpgrid(uu,dya=dya)
		vv = interpgrid(vv,dya=dya)

	
	if interp:
	
		uu2 = uu
		vv2 = vv
		
	else:
	
		uu2 = OrderedDict()
		vv2 = OrderedDict()

		for start,end in zip(dist[0:-1],dist[1::]):
			
			uu2.update({np.mean([start,end]).astype('string'): uu.iloc[:,(uu.columns.values>start)&(uu.columns.values<end)].mean(axis=1).values})
			vv2.update({np.mean([start,end]).astype('string'): vv.iloc[:,(vv.columns.values>start)&(vv.columns.values<end)].mean(axis=1).values})


		uu2 = pd.DataFrame.from_dict(uu2)
		vv2 = pd.DataFrame.from_dict(vv2)


	ang = np.mean(ang)
#	cosang = np.cos(ang)
	print ang
	if ang>0:
		vel_obs = (vv2*np.cos(ang))-(uu2*np.sin(ang))
		vel_brute = (vv*np.cos(ang))-(uu*np.sin(ang))
	elif ang<0:
		vel_obs = (vv2*np.cos(ang))-(uu2*np.sin(ang))
		vel_brute = (vv*np.cos(ang))-(uu*np.sin(ang))
	else:
		print 'Problema no Ângulo médio da radial'

	plt.figure()
	plt.contourf(distv,
	-1*np.arange(0,601),
	vel_brute,100,
	cmap='RdBu_r',vmin=-1.5,vmax=1.5)
	plt.colorbar()
	c2 = plt.contour(distv,
	-1*np.arange(0,601),
	vel_brute,colors='k',levels=np.arange(-1.5,1.5,0.1))
	clabels = plt.clabel(c2, fontsize=10,fmt='%1.2f',
	colors='k',use_clabeltext=True)
	plt.plot(dist,np.repeat(0,dist.shape),
	marker='v',markersize=20,color='k',
	linewidth=0)

	try:
		vel_obs = pd.DataFrame(vel_obs.values,index=np.arange(0,dya+1),columns=vv2.columns)	
	except:
		vel_obs = pd.DataFrame(vel_obs.values,index=np.hstack(xy['zc']).astype('int'),columns=vv2.columns)


	
	return latv,lonv,vel_obs,yy,xx,vel_brute,uu,vv


	
def woa_sec(secref,path = '/home/iury/Copy/TCC/Rotinas/dados'):
	'''
	This function extract from woa climatology the nearest
	profiles from some CTD section as pandas.Panel
	
	Uses the function interpgrid from seacalc
	'''
	
	#load the data woa.mat 
	woa = sio.loadmat(os.path.join(path,'woa.mat'))
	#define variables of climatology
	tclim = woa['tclim']
	sclim = woa['sclim']
	
	#define coordinates from section
	lonsta = secref.minor_xs('lon').mean()
	latsta = secref.minor_xs('lat').mean()
	
	#find the equivalent coordinates from woa data
	lons = np.round(lonsta).values+360
	lats = np.round(latsta).values+90
	
	#get the temperature and salinity of these points
	tt = tclim[lats.tolist(),lons.tolist(),:]
	ss = sclim[lats.tolist(),lons.tolist(),:]

	#pre-define the depth levels of woa data
	woa_dpt = np.array([0,10,20,30,50,75,100,125,150,200,250,300,400,500,
	         600,700,800,900,1000,1100,1200,1300,1400,1500,
	         1750,2000,2500,3000,3500,4000,4500,5000,5500])
	
	#calculate the distance vector from real data points
	dist = np.hstack((0,np.cumsum(sw.dist(latsta,lonsta)[0])))
	
	#create pandas.DataFrame for each property
	tt = pd.DataFrame(tt.T,index=woa_dpt,columns=dist)
	ss = pd.DataFrame(ss.T,index=woa_dpt,columns=dist)
	
	#interpolate the data to 1m step of depth
	tt = interpgrid(tt,dya=1000)
	ss = interpgrid(ss,dya=1000)
	
	return tt,ss


	
def etopo_sec(lat,lon,interp=False,path = '/home/iury/Copy/TCC/Rotinas'):
	'''
	This function returns the section data from etopo bathymetry to plot
	it on section contours.
	path = path to dir of .mat data from etopo
	etopo data must be in .mat format as an struct  'bat' with x,y and z
	'''
	
	#reading etopo data
	etopo = sio.loadmat(os.path.join(path,'etopo_tcc.mat'))

        #define the variables. [0][0] is necessary because when transform
        #the .mat data to python it puts a lot of brackets
	batx = etopo['bat']['x'][0][0].astype('float')
	baty = etopo['bat']['y'][0][0].astype('float') #convert to float
	Batx,Baty = np.meshgrid(batx,baty) # create a grid
	batz = etopo['bat']['z'][0][0].astype('float')
	

        # extracting the bathymetry section based on latitudes and longitudes
	batsect = []
	for lat1,lon1 in zip(lat,lon):
	        #find the minimum difference between data longitudes
	        #and etopo longitudes
		xpos = np.argmin(np.abs(batx-lon1)) # abs function return the module
		ypos = np.argmin(np.abs(baty-lat1)) 
		#append function add an element to the end of a list
		batsect.append(batz[ypos,xpos])
	
	#the sw.dist returns distance and angles, [0] take only the distances
	dist = np.hstack((0,np.cumsum(sw.dist(lat,lon)[0])))# hstack concatenate all values
	
	# if interp is true the function interpolates the data to smooth bathymetry
	if interp:
		#new distance vector from 0 to maxmimum with 0.1 of step
		new_dist = np.arange(0,dist.max(),0.1)
                #f is the interpolation function from etopo data
		f = scint.interp1d(dist,batsect,kind='cubic')
		#uses the f function to get the new bathymetry data
		new_bat = f(new_dist)
		#convert it to numpy array
		batsect = np.array(batsect)
		#create a vector of repeated value less than minimum of bathymetry
		#this is used to create the polygon of bathymetry on image
		d = np.tile(new_bat.min()-10,new_bat.shape)
		
		return new_dist,new_bat
	else:
		return dist,batsect
	

def extrap_section(lon,lat,section):
	'''
	This function extrapolate the properties from CTD 
	pandas.Panel section using the gradient and the functions
	and extrap_gradient
	
	(numpy.array,numpy.array,pandas.Panel) => pandas.Panel
	'''

	#counter of atributes
	cnt = 0
	#for each atribute
	for atribute in section.minor_axis:
		#select the pandas.DataFrame of atribute
		atrb = section.minor_xs(atribute)
		#extrapolate the data using the gradient
		atrb = extrap_gradient(atrb,coords=True,lat=lat,lon=lon)
		#substitute the values of atribute in section
		section.iloc[:,:,cnt] = atrb
		#add 1 to counter
		cnt += 1
	
	return section	


