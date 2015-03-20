# -*- coding: utf-8 -*-
#Most of these functions is used to extract data from raw files
#and organize this data to some better format
#Made by Iury Sousa - São Paulo/Brazil

#Some of these functions were based or inspired by other authors
#as Filipe Fernandes and his scripts published on GitHub or his 
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


import gsw
import seawater as sw


#################################################################


def return_section(directory):
	'''
	This function reads all the files from directory
	that has no extension as a pandas pickle, sorted
	by name and return a panel section and its lat 
	and lon.
	'''

	lon,lat = [],[] # latitude and longitude starts empty

	# this is the list of cnv files from directory
	# sorted by name
	stas = np.sort(glob(os.path.join(directory,'*.cnv')))
	# define section as an empty ordered dictionary
	section = OrderedDict()
	# starts the loop
	for pf in stas:
	# read the file with the same name and no extension
	# because pandas pickle has no extensio
		data = pd.read_pickle(pf.split('.')[0])
	# define the name of the file for the dict
	# again, the name has no extension (use split)
		name = os.path.basename(pf).split('.')[0]
	# give to section an update, so this dict
	# will have a variable with the same file name
		section.update({name: data})
	# append to lon and lat
		lon.append(data.lon.mean())
	# the mean is because in data the lat and lon is
	# a column with repeated lat and lon
		lat.append(data.lat.mean())
		
	if np.shape(section[name])[0]<1800:
		section = OrderedDict(section.items()[::-1])
		lat.reverse()
		lon.reverse()
	# transform from dict to panel
	section = pd.Panel.fromDict(section)
	
	section = section.fillna(method='backfill') 
	
	#if the last station have mone nans than the first
	if section.minor_xs('t090').isnull().sum()[[1,-1]].diff()[-1]>0:
	    #invert the items (stations) dimension
	    section = section.iloc[::-1,:,:]
	
	# return the lat, lon and section
	return lat,lon,section
	




def extsect(radpath,isopicnal=32.15,
	prof=False,refer=1000,sup=False):
	'''
	This function reads adn create a geostrophic
	velocity section reffered on some isopicnal sigma_1
	'''

	lat,lon,dist,pref,section = [],[],[],[],[]
	# reads the files and return section lats and lons
	# using return_section function
	
	section = pd.read_pickle(os.path.join(radpath,os.path.split(radpath)[1]))
	
	lat = section.minor_xs('lat').mean().values
	lon = section.minor_xs('lon').mean().values

	
	# defines the radial name and print it
	radname = 'r'+os.path.split(radpath)[1]
	print radname

	# define distance vector for the section
	dist = np.cumsum(np.append(0, sw.dist(lat, lon,
	units='km')[0]))
	# define distance vector the de middle point
	# in section, used for geostrophic velocity
	distv = np.cumsum(np.append(0, sw.dist(lat, lon,
	units='km')[0][0:-1]))+sw.dist(lat, lon, units='km')[0]/2

	# create a figure
	plt.figure()
	# make filled contour plot of sigma_1
	plt.contourf(dist,-section.major_axis,
	section.minor_xs('psigma1'),100, cmap='jet')

	# make a contour plot of sigma_1 on the selected
	# isopicnal with handle c1
	c1 = plt.contour(dist,-section.major_axis,
              section.minor_xs('psigma1'),
	      np.array([isopicnal]))
	# find the coordinates of the contour plot
	# vstack is used because u may have many arrays
	# in allsegs because the contour may be broken
	coords = np.vstack(c1.allsegs[0])

	# find the distance and it depth (as integer)
	# for the isopicnal
	x = coords[:,0]
	y = np.round(coords[:,1]).astype('int')

#	plt.close()



	# find the the depth of isopicnal for each point in dist
	for i in dist:
	# if the profile has no data in this depth
	    if np.array([y[i==x]]).size == 0:
		    # uses 1600 as reference (this value is used to convert profile in nan)
		    # because the 32.15 sigma1 is in less depth from this
		    # so the velocity in 1600 will be nan
	            pref.append(1600)
	# if the profile has data in this depth
	    else:
	# uses the found reference depth
		    pref.append(np.mean(-y[i==x]).astype('int'))
	
	# the reference depth is for distance for stations
	# it must be converted for the distance between the middle
	# points, in other words, for velocity section	
	pref = np.round(pref[0:-1]+np.diff(pref)/2)
	# geopotential section
	gpan = section.minor_xs('gpan')
	# if we have the shallowest profile in the end, it must be
	# inverted 
	if pref[-1]==1600:
			print radname+' está invertida!'
	# invert lats and lons
			lat = np.flipud(lat)
			lon = np.flipud(lon)
	# define the new dist and distv
			dist = np.cumsum(np.append(0, sw.dist(lat, lon,
			units='km')[0]))
			distv = np.cumsum(np.append(0, sw.dist(lat, lon,
			units='km')[0][0:-1]))+sw.dist(lat, lon, units='km')[0]/2
	# invert pref and gpan
			pref = np.flipud(pref)
			gpan =  np.fliplr(gpan)

	# calculates the geostrophic velocity based on surface
	velsup = sw.gvel(gpan,lat,lon)
	velsup = np.asmatrix(velsup)


	if prof:
		# if u want 1000 meters referenced geostrophic velocity
		vel_iso = np.subtract(velsup,velsup[refer,:])	
		
	elif sup:
	
		vel_iso = velsup
		
	
	else:
		# reference the velocity to the isopicnal
		vel_iso = np.subtract(velsup,velsup[pref.astype('int'),np.arange(0,velsup.shape[1])])
	


	# interpolate the velocity because some depths we have no data 
	if interp:
		vel_iso = pd.DataFrame(vel_iso).interpolate(method='linear', axis=0).values

	vel_iso = pd.DataFrame(vel_iso,index=section.major_axis,columns=distv)
	# print the end of the extrated radial
	print '------------------'
	# return the x, y, z, lats and lons
	return vel_iso,lat,lon,section
	
	



def extvel(section,pathadcp,
	refer=1000,interp=False,
	adcpv=False,fmed=True):

	pref = []

	lat = section.minor_xs('lat').mean().values
	lon = section.minor_xs('lon').mean().values
	
#	lata = np.linspace(lat[0],lat[-1],20)
#	lona = np.linspace(lon[0],lon[-1],20)
	latmid = lat[:-1] + np.diff(lat)/2
	lonmid = lon[:-1] + np.diff(lon)/2

	if interp:	
		lata,lona = [],[]
		for i in np.arange(0,lat.size):
		
			try:
				lata.append(lat[i])
				lona.append(lon[i])
				lata.append(latmid[i])
				lona.append(lonmid[i])
			except:
				print('')
	else:
		lata,lona = lat,lon			
			

	lat_adcp,lon_adcp,adcp,lat_brute,lon_brute,vel_brute,vv,uu = adcpvel(lata,lona,pathadcp,dya=600,interp=adcpv)




	
	# defines the radial name and print it
#	radname = 'r'+os.path.split(radpath)[1]
#	print radname


	dista = np.cumsum(np.append(0,
	sw.dist(lata, lona)[0]))
	
	distva = dista[:-1] + np.diff(dista)/2

	mask = ~np.all(np.isnan(adcp),axis=0)
	distva = distva[mask]
	adcp = adcp.loc[:,mask]

	# define distance vector for the section
	dist = np.cumsum(np.append(0,
	sw.dist(lat, lon)[0]))
	
	# define distance vector the de middle point
	# in section, used for geostrophic velocity
	distv = np.cumsum(np.append(0, sw.dist(lat, lon,
	units='km')[0][0:-1]))+sw.dist(lat, lon, units='km')[0]/2


#	
#	# the reference depth is for distance for stations
#	# it must be converted for the distance between the middle
#	# points, in other words, for velocity section	
#	pref = np.round(pref[0:-1]+np.diff(pref)/2)

	# geopotential section
	gpan = section.minor_xs('gpan')
	gpan = pd.DataFrame(section.minor_xs('gpan').values,columns=dist)
	
	if interp:
		gpan = innerinterp(gpan,dista)

	adcp = pd.DataFrame(adcp.values,columns=distva)
	adcpmdr = adcp.copy()
	adcpmdr.loc[adcpmdr.index.values<100] = np.nan
	
	
	# calculates the geostrophic velocity based on surface
	velsup = sw.gvel(gpan,lata,lona)
	velsup = np.asmatrix(velsup)

	velsup = velsup[:,mask]

#	# 1000 meters referenced geostrophic velocity
#	vel_prof = np.subtract(velsup,velsup[refer,:])	
#	

##	converting to dataframes
	velsup = pd.DataFrame(velsup,columns=distva)
	latva = lata[:-1] + np.diff(lata)/2
	if fmed:
		f = sw.f(lata[mask])
		velsup = velsup*(f/np.mean(f))
#	vel_prof = pd.DataFrame(vel_prof,index=section.major_axis,columns=distv)
#	vel_iso = pd.DataFrame(vel_iso,index=section.major_axis,columns=distv)

#	referenced to adcp data geostrophic velocity
	vel_mdr,refmdr = mdr(velsup,adcpmdr)
#	vel_mdr,refmdr = mdr(velsup,adcp)
	vel_mdr = pd.DataFrame(vel_mdr,columns=distva)

#	referenced to adcp data at different depths geostrophic velocity		
	vel_fixmdr = mdr_fixed(velsup,adcp,150)
	vel_fixmdr = pd.DataFrame(vel_fixmdr,columns=distva)



	prof = adcp.index.values
	LATV,y = np.meshgrid(lat_adcp,prof)
	LONV,y = np.meshgrid(lon_adcp,prof)
	


	dens0 = pd.DataFrame(section.minor_xs('psigma0').values,columns=dist)
	dens1 = pd.DataFrame(section.minor_xs('psigma1').values,columns=dist)
	
	t090 = pd.DataFrame(section.minor_xs('t090').values,columns=dist)
	sp = pd.DataFrame(section.minor_xs('sp').values,columns=dist)

	if interp:
		dens0 = innerinterp(dens0,distva)
		t090 = innerinterp(t090,distva)
		dens1 = innerinterp(dens1,distva)
		sp = innerinterp(sp,distva)
	else:
		dens0 = dens0.iloc[:,:-1]+np.diff(dens0.values)/2
		dens1 = dens1.iloc[:,:-1]+np.diff(dens1.values)/2
		t090 = t090.iloc[:,:-1]+np.diff(t090.values)/2
		sp = sp.iloc[:,:-1]+np.diff(sp.values)/2
		
		dens0 = pd.DataFrame(dens0.loc[:,mask].values,columns=distva)
		dens1 = pd.DataFrame(dens1.loc[:,mask].values,columns=distva)
		t090 = pd.DataFrame(t090.loc[:,mask].values,columns=distva)
		sp = pd.DataFrame(sp.loc[:,mask].values,columns=distva)
		
		

	prof = dens0.index.values
	
	LATV,y = np.meshgrid(lat_adcp[mask],prof)
	LONV,y = np.meshgrid(lon_adcp[mask],prof)
	
	print t090.shape
	print sp.shape
	print y.shape

	dens2 = sw.pden(t090.values,sp.values,y,pr=2000)
	dens4 = sw.pden(t090.values,sp.values,y,pr=4000)

	dens2 = pd.DataFrame(dens2,columns=distva)
	dens4 = pd.DataFrame(dens4,columns=distva)

	LATV = pd.DataFrame(LATV,index=velsup.index.values,columns=distva)
	LONV = pd.DataFrame(LONV,index=velsup.index.values,columns=distva)

	
	secvel = OrderedDict()
	secvel.update({'dens0':dens0})
	secvel.update({'dens1':dens1})
	secvel.update({'dens2':dens2})
	secvel.update({'dens4':dens4})
	secvel.update({'velsup':velsup})
	secvel.update({'velmdr':vel_mdr})
	secvel.update({'velfixmdr':vel_fixmdr})
#	secvel.update({'velprof':vel_prof})
	secvel.update({'adcp':adcp})
	secvel.update({'lat':LATV})
	secvel.update({'lon':LONV})
		
	secvel = pd.Panel.from_dict(secvel)
	
	return secvel,refmdr,vel_brute



def extvel2(section,pathadcp,isopicnal=32.15,refer=1000,interpadcp=False):

	pref = []

	lat = section.minor_xs('lat').mean().values
	lon = section.minor_xs('lon').mean().values

	lat_adcp,lon_adcp,adcp,lat_brute,lon_brute,vel_brute,vv,uu = adcpvel(lat,lon,pathadcp,dya=600,interp=interpadcp)
	
	# defines the radial name and print it
#	radname = 'r'+os.path.split(radpath)[1]
#	print radname

	# define distance vector for the section
	dist = np.cumsum(np.append(0,
	sw.dist(lat, lon)[0]))
	
	# define distance vector the de middle point
	# in section, used for geostrophic velocity
	distv = np.cumsum(np.append(0, sw.dist(lat, lon,
	units='km')[0][0:-1]))+sw.dist(lat, lon, units='km')[0]/2


	# create a figure
	plt.figure()
	# make filled contour plot of sigma_1
#	plt.contourf(dist,-section.major_axis,
#	section.minor_xs('psigma1'),100, cmap='jet')

	# make a contour plot of sigma_1 on the selected
	# isopicnal with handle c1
	c1 = plt.contour(dist,-section.major_axis,
              section.minor_xs('psigma1'),
	      np.array([isopicnal]))
	# find the coordinates of the contour plot
	# vstack is used because u may have many arrays
	# in allsegs because the contour may be broken
	coords = np.vstack(c1.allsegs[0])

	# find the distance and it depth (as integer)
	# for the isopicnal
	x = coords[:,0]
	y = np.round(coords[:,1]).astype('int')

	plt.close()



	# find the the depth of isopicnal for each point in dist
	for i in dist:
	
#	TALVEZ ISSO NÃO SEJA MAIS NECESSÁRIO JÁ QUE AS SEÇÕES FORAM EXTRAPOLADAS
	# if the profile has no data in this depth
	    if np.array([y[i==x]]).size == 0:
		    # uses 1600 as reference (this value is used to convert profile in nan)
		    # because the 32.15 sigma1 is in less depth from this
		    # so the velocity in 1600 will be nan
	            pref.append(1600)
	# if the profile has data in this depth
	    else:
	# uses the found reference depth
		    pref.append(np.mean(-y[i==x]).astype('int'))

	
	# the reference depth is for distance for stations
	# it must be converted for the distance between the middle
	# points, in other words, for velocity section	
	pref = np.round(pref[0:-1]+np.diff(pref)/2)

	# geopotential section
	gpan = section.minor_xs('gpan')

	# calculates the geostrophic velocity based on surface
	velsup = sw.gvel(gpan,lat,lon)
	velsup = np.asmatrix(velsup)



	# 1000 meters referenced geostrophic velocity
	vel_prof = np.subtract(velsup,velsup[refer,:])	
	
		
	# reference the velocity to the isopicnal
	vel_iso = np.subtract(velsup,velsup[pref.astype('int'),np.arange(0,velsup.shape[1])])


#	converting to dataframes
	velsup = pd.DataFrame(velsup,index=section.major_axis,columns=distv)
	vel_prof = pd.DataFrame(vel_prof,index=section.major_axis,columns=distv)
	vel_iso = pd.DataFrame(vel_iso,index=section.major_axis,columns=distv)

#	referenced to adcp data geostrophic velocity
	vel_mdr,refmdr = mdr(velsup,adcp)
	vel_mdr = pd.DataFrame(vel_mdr,index=section.major_axis,columns=distv)

#	referenced to adcp data at different depths geostrophic velocity		
	vel_fixmdr = mdr_fixed(velsup,adcp,150)
	vel_fixmdr = pd.DataFrame(vel_fixmdr,index=section.major_axis,columns=distv)

	adcp = pd.DataFrame(adcp.values,columns=velsup.columns.values.astype('float'))

#	uncomment if u want that she surface nans be extrapolated
#	adcp.iloc[:] = np.flipud(adcp)
#	adcp = adcp.interpolate(method='linear',axis=0)
#	adcp.iloc[:] = np.flipud(adcp)


	prof = adcp.index.values
	LATV,y = np.meshgrid(lat_adcp,prof)
	LONV,y = np.meshgrid(lon_adcp,prof)
	

	dens0 = section.minor_xs('psigma0')
	dens1 = section.minor_xs('psigma1')
	
	t090 = section.minor_xs('t090')
	sp = section.minor_xs('sp')
	
	
	dens0 = dens0.iloc[:,:-1]+np.diff(dens0.values)/2
	dens1 = dens1.iloc[:,:-1]+np.diff(dens1.values)/2
	t090 = t090.iloc[:,:-1]+np.diff(t090.values)/2
	sp = sp.iloc[:,:-1]+np.diff(sp.values)/2


	dens0 = pd.DataFrame(dens0.values,index=section.major_axis,columns=distv)
	dens1 = pd.DataFrame(dens1.values,index=section.major_axis,columns=distv)

	prof = dens0.index.values
	LATV,y = np.meshgrid(lat_adcp,prof)
	LONV,y = np.meshgrid(lon_adcp,prof)

	dens2 = sw.pden(t090.values,sp.values,y,pr=2000)
	dens4 = sw.pden(t090.values,sp.values,y,pr=4000)

	dens2 = pd.DataFrame(dens2,index=section.major_axis,columns=distv)
	dens4 = pd.DataFrame(dens4,index=section.major_axis,columns=distv)

	LATV = pd.DataFrame(LATV,index=velsup.index.values,columns=distv)
	LONV = pd.DataFrame(LONV,index=velsup.index.values,columns=distv)

	
	secvel = OrderedDict()
	secvel.update({'dens0':dens0})
	secvel.update({'dens1':dens1})
	secvel.update({'dens2':dens2})
	secvel.update({'dens4':dens4})
	secvel.update({'velsup':velsup})
	secvel.update({'velmdr':vel_mdr})
	secvel.update({'velfixmdr':vel_fixmdr})
	secvel.update({'velprof':vel_prof})
	secvel.update({'veliso':vel_iso})
	secvel.update({'adcp':adcp})
	secvel.update({'lat':LATV})
	secvel.update({'lon':LONV})
		
	secvel = pd.Panel.from_dict(secvel)
	
	return secvel,refmdr,vel_brute