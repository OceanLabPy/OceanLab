# -*- coding: utf-8 -*-
#Most of these functions is used to plot the data extracted,
#calculated and organized by the functions from seacalc and seaext
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
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
from matplotlib import gridspec
from matplotlib.transforms import Bbox
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm

#import fnmatch


import gsw
import seawater as sw
import matplotlib.patches
import matplotlib



def dfcontourf(ax,df,levels=None,vmin=None,vmax=None,cmap='jet',ocean_pressure=True):
        '''
        This function does filled contour plot using index and columns as axes from
        DataFrame pandas object.
        To be useful the df must be with axis values as index and columns.
        '''
        #define x as the columns values
	x = df.columns.values
	#if the index represent pressure in ocean
	if ocean_pressure:
          y = df.index.values*-1. # convert it to negative float
        else:
          y = df.index.values*1. # convert it to float
	
	#create grid with coordinates
	X,Y = np.meshgrid(x,y)
	#define Z as the DataFrame values
	Z = df.values

        #define vmin,vmax,and levels if it hasn't defined already
	if vmin==None:
		vmin=np.nanmin(Z)
	if vmax==None:
		vmax=np.nanmax(Z)
	if levels==None:
		levels=np.arange(vmin,vmax,0.1)
	
	#print the extreme of values to colormap
	print vmax
	print vmin

        #create the contour with h as a handle or object
	h = ax.contourf(X,Y,Z,levels,
		vmin=vmin,vmax=vmax,
		cmap=cmap)
        #return the contourf object
	return h
	
def dfcontour(ax,df,levels=None,color='k',ocean_pressure=True):
        '''
        This function does filled contour plot using index and columns as axes from
        DataFrame pandas object.
        To be useful the df must be with axis values as index and columns.
        '''
        #define x as the columns values
	x = df.columns.values
        #if the index represent pressure in ocean
        if ocean_pressure:
          y = df.index.values*-1. # convert it to negative float
        else:
          y = df.index.values*1. # convert it to float
	
        #create grid with coordinates
        X,Y = np.meshgrid(x,y)
        #define Z as the DataFrame values
        Z = df.values

        #define vmin,vmax,and levels if it hasn't defined already
        if vmin==None:
                vmin=np.nanmin(Z)
        if vmax==None:
                vmax=np.nanmax(Z)
	if levels==None:
		levels=np.arange(vmin,vmax,0.1)
		
	print vmax
	print vmin

        #create the contour with h as a handle or object
	h = ax.contour(X,Y,Z,levels,
		colors=color)
        #create the contour label of values as a float with .1 of precision
	ax.clabel(h, fontsize=10, fmt='%.1f')
	#return the object h
	return h
	
	
def velcomp(path,left=0.8,mmax = 1.5):
	data  = pd.read_pickle(os.path.join(path,'data'))
	velmdr = data['secvel'].velmdr
	adcp = data['secvel'].adcp


	secvel = data['secvel']
	lon = secvel.lon.mean().values
	lat = secvel.lat.mean().values


	a,bathy = etopo_sec(lat,lon)
	bathy = np.array(bathy)

	d = np.tile(bathy.min()-10,bathy.shape)

	fig,(ax1,ax2) = plt.subplots(2,1,figsize=(10,10))

	levs = np.arange(-mmax,mmax+0.1,0.1)
	distmdr = velmdr.columns.values
	distadcp = adcp.columns.values

	mask = ~np.isnan(adcp.loc[150].values)
	c1 = dfcontourf(ax1,velmdr,levs,
		vmin=-mmax,vmax=mmax,cmap='RdBu_r')
	c11 = dfcontour(ax1,velmdr,levs)
	ax1.set_xlim([adcp.columns.values[mask][0],adcp.columns.values[mask][-1]])
	ax1.set_ylim([-400,0])
	ax1.plot(distmdr,-data['refmdr'],'r',linewidth=2)
	ax1.plot(distmdr,np.tile(-1,distmdr.shape),
		color='k',marker = 'v',
		linewidth=0,markersize=20)
	ax1.grid('on')
	ax1.set_xticklabels([])

	ax1.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)

	c2 = dfcontourf(ax2,adcp,levs,
		vmin=-mmax,vmax=mmax,cmap='RdBu_r')
	c22 = dfcontour(ax2,adcp,levs)
	ax2.set_xlim([adcp.columns.values[mask][0],adcp.columns.values[mask][-1]])
	ax2.plot(distadcp,np.tile(-1,distadcp.shape),
		color='k',marker = 'v',
		linewidth=0,markersize=20)
	ax2.set_ylim([-400,0])
	ax2.set_xlabel(u'Distância',fontweight='bold')
	ax2.grid('on')

	ax2.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)
	
	position=fig.add_axes([left,0.1,0.02,0.8])
	cbar = fig.colorbar(c2,ticks=np.arange(-mmax,mmax+0.1,0.2),cax=position)

	return fig,ax1,ax2,cbar



def interpvel(lons,lats,xx,yy,u1,v1):
	points = zip(lons,lats)
	newpoints = zip(xx.ravel(),yy.ravel())
	
	u = scint.griddata(points,u1,newpoints)
	v = scint.griddata(points,v1,newpoints)
	
	u = u.reshape(xx.shape)
	v = v.reshape(xx.shape)
	
	return u,v
	
def streams(ax,fig,xx,yy,u,v,vmax=1.5,left=0.8):
	'''
	This function create streamplot for curvilinear
	evenly spaced grid for some given axis and figure.
	
	streams(ax,fig,xx,yy,u,v) -> streamplot
	'''
	x = np.linspace(xx.min(), xx.max(), 50)
	y = np.linspace(yy.min(), yy.max(), 50)

	speed = np.sqrt((u**2)+(v**2))

	xi, yi = np.meshgrid(x,y)

	#then, interpolate your data onto this grid:
	
	px = xx.flatten()
	py = yy.flatten()
	pu = u.flatten()
	pv = v.flatten()
	pspeed = speed.flatten()
	
	gu = griddata(zip(px,py), pu, (xi,yi))
	gv = griddata(zip(px,py), pv, (xi,yi))
	gspeed = griddata(zip(px,py), pspeed, (xi,yi))
	
	lw = 4*gspeed/vmax
	#now, you can use x, y, gu, gv and gspeed in streamplot:
	
	
	xx,yy = ax(xx,yy)
	xi,yi = ax(xi,yi)
	
#	ax.contour(xx,yy,speed, colors='k', alpha=0.4)
	ax.plot(xx,yy,'-k',alpha=0.1)
	ax.plot(xx.T,yy.T,'-k',alpha=0.1)
#	ax.plot(xi,yi,'-b',alpha=0.1)
#	ax.plot(xi.T,yi.T,'-b',alpha=0.1)


	### Create a list of RGB tuples
	colors = [(0.847, 0.057, 0.057),
	         (0, 0.592, 0),
	         (0, 0.559, 0.559),
	         (0.316, 0.316, 0.991),
	         (0.718, 0, 0.718)] # This example uses the 8-bit
	colors = np.flipud(colors).tolist()
	### Call the function make_cmap which returns your colormap
	my_cmap = make_cmap(colors, bit=False)


	c = ax.streamplot(xi,yi,gu,gv, density=2,
		linewidth=lw, color=gspeed, cmap=my_cmap)
	
		
	ax.streamplot(xi,yi,gu,gv, density=1,
		linewidth=lw, color='k')

	c.arrows.set_clim([0,vmax])
	c.lines.set_clim([0,vmax])		
	position=fig.add_axes([left,0.1,0.02,0.8])
	cbar = fig.colorbar(c.lines,cax=position)
	
	plt.draw()



	# iterate through the children of ax
	for art in ax.ax.get_children():
	    # we are only interested in FancyArrowPatches
	    if not isinstance(art, matplotlib.patches.FancyArrowPatch):
	        continue
	    # remove the edge, fill with black
	    art.set_edgecolor([0, 0, 0, 1])
	    art.set_facecolor([0, 0, 0, 1])
#	    # make it bigger
#	    art.set_mutation_scale(20)
	    # move the arrow head to the front
	    art.set_zorder(10)
#	cbar.set_clim([0,1.5])

#	c.arrows.set_cmap(None)
#	c.arrows.set_color('k')

#	c = ax.streamplot(xi,yi,gu,gv, density=3,
#		linewidth=lw, color='k')
	return c
	
	
#################################################################


def make_map(llcrnrlon=-51, urcrnrlon=-30.2, llcrnrlat=-15,
		urcrnrlat=7.1,projection='merc', resolution='i',
		figsize=(6, 6), inset=True,axe = None,steplat=2,
		steplon=2,inloc=1,contcolor='0.85'):
    
    '''
    This function creates a basemap map easily with option
    of inset axe with geolocation.
    
    It's all based on make_map from Filipe Fernandes, but with some
    changes
    '''



    if axe==None:
            fig, ax = plt.subplots(figsize=figsize)
	    m = Basemap(llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                projection=projection, resolution=resolution)
    else:
	    m = Basemap(llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                projection=projection, resolution=resolution,ax=axe)
	    ax = axe
            print('Map Created!')




    m.drawstates(zorder=22)
    m.drawcoastlines(zorder=21)
    m.fillcontinents(color=contcolor,zorder=20)
    meridians = np.arange(np.floor(llcrnrlon), 
		np.ceil(urcrnrlon) + 2, steplon)
    parallels = np.arange(np.floor(llcrnrlat), 
		np.ceil(urcrnrlat) + 2, steplat)
    m.drawparallels(parallels, linewidth=0.4, labels=[1, 0, 0, 0],zorder=23)
    m.drawmeridians(meridians, linewidth=0.4, labels=[0, 0, 0, 1],zorder=24)
    m.llcrnrlon = llcrnrlon
    m.urcrnrlon = urcrnrlon
    m.llcrnrlat = llcrnrlat
    m.urcrnrlat = urcrnrlat
    m.ax = ax

    if inset:
        axin = inset_axes(m.ax, width="30%", height="30%", loc=inloc)
        # Global inset map.
        inmap = Basemap(projection='ortho', 
	lon_0=np.mean(-42), lat_0=np.mean(-4),
        ax=axin, anchor='NE')
        inmap.drawcountries(color='white')
        inmap.fillcontinents(color='gray')
        bx, by = inmap(m.boundarylons, m.boundarylats)
        xy = list(zip(bx, by))
        mapboundary = Polygon(xy, edgecolor='r', linewidth=1, fill=False)
        inmap.ax.add_patch(mapboundary)
    if axe==None:
	    return fig, m
    else:
	    return m


def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
    

def tscomp(secref):
	'''
	This function creates T-S plots for CTD stations
	and referred climatology from WOA (NOAA).
	'''
	path = '/home/iury/Copy/TCC/dados/CLIMAT/WOA'

	woa = sio.loadmat(os.path.join(path,'woa.mat'))

	tclim = woa['tclim']
	sclim = woa['sclim']
	
	lons = np.round(secref.minor_xs('lon').mean()).values+360
	lats = np.round(secref.minor_xs('lat').mean()).values+90
	
	tt = tclim[lats.tolist(),lons.tolist(),:].ravel()
	ss = sclim[lats.tolist(),lons.tolist(),:].ravel()
	
	fig,ax = plt.subplots()
	

	colors = cm.rainbow(np.linspace(0, 1, len(secref.items)))
	
	ax.plot(ss,tt,label='WOA',marker='o',color='k',linewidth=0,markersize=5)
	for i,c in zip(secref.items,colors):
		ax.scatter(secref.minor_xs('sp')[i],secref.minor_xs('pt')[i],
		label=i,alpha=0.8,edgecolor='None',color=c,s=7)
		
	ax.legend(fontsize=10,loc=4)
	
	clev = [20,25.6,26.9,27.38,27.53,28]
	
	sal,temp = np.meshgrid(np.arange(33,39),np.arange(0,31))
	dens = sw.dens0(sal,temp)-1000
	dens1 = sw.dens(sal,temp,1000)-1000
	dens2 = sw.dens(sal,temp,2000)-1000
#	dens4 = sw.dens(sal,temp,4000)-1000

	
#	c1 = ax.contourf(sal,temp,dens,levels=clev,alpha=0.5,
#		colors=('#ff0000', '#ff9900', '#00c4ff', '#4fda66', '#0000ff'))
	
	c2 = ax.contour(sal,temp,dens,levels=[24.5,26.8],colors='k')	
	plt.clabel(c2, fontsize=10)
	c2 = ax.contour(sal,temp,dens1,levels=[32.15],colors='k')	
	plt.clabel(c2, fontsize=10)
	c2 = ax.contour(sal,temp,dens2,levels=[37],colors='k')	
	plt.clabel(c2, fontsize=10)
	
	plt.xlim([sal.min().min(),sal.max().max()])
	plt.ylim([temp.min().min(),temp.max().max()])
	plt.xlabel('Salinidade')
	plt.ylabel('Temperatura')
#	c2 = ax.contour(sal,temp,dens4,levels=[45.83,45.9],colors='k')	
#	plt.clabel(c2, fontsize=10)
	
			
#	position=fig.add_axes([0.11,0.06,0.8,0.02])
#	fig.colorbar(c1,orientation='horizontal',cax = position)
#	
#	position.text(x=0.1,y=-1.5,ha='center',s='AT')
#	position.text(x=0.3,y=-1.5,ha='center',s='ACAS')
#	position.text(x=0.5,y=-1.5,ha='center',s='AIA')
#	position.text(x=0.7,y=-1.5,ha='center',s='ACS')
#	position.text(x=0.9,y=-1.5,ha='center',s='APAN')
	
	pos = np.array(ax.get_position().bounds)
	pos[1] += pos[1]*0.35
	ax.set_position(pos)
	ax.set_title('Diagrama TS')
	
	return fig,ax
	
	
	

def plotvsec(secvel,lim=-500):
	exec('vsec = secvel.vel'+velname)
	dens0 = secvel.dens0
	
	X,Y = np.meshgrid(vsec.columns.values,-vsec.index.values)
	vlevs = np.arange(-1.5,1.5,0.1)

	fig,ax = plt.subplots()
	ax.contourf(X,Y,vsec,100,cmap='RdBu_r',vmin=-1.5,vmax=1.5)
	c1 = ax.contour(X,Y,vsec,levels=vlevs,colors='k')
	
	ax.clabel(c1,fmt='%1.1f', fontsize=10,inline_spacing=0)

	X,Y = np.meshgrid(dens0.columns.values.astype('float'),-dens0.index.values)
	dlevs = np.array([24.5,26.5,26.75,27])
	
	c2 = ax.contour(X,Y,dens0,levels=dlevs,colors='k',
	linewidths=2,linestyles='dashed')
	
	clabels = ax.clabel(c2, fontsize=10,fmt='%1.2f',
	colors='k',use_clabeltext=True)
	
	plt.ylim([lim,0])
	

	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))
		
		
def plotadcp(secvel,lim=-400,radcomplicada=False):
	adcp = secvel.adcp
	dens0 = secvel.dens0
	
	lon = secvel.lon.mean().values
	lat = secvel.lat.mean().values

	if radcomplicada:
		lat[3] = -2.0933934592636545
		lon[3] = -38.442823159524877

	a,bathy = etopo_sec(lat,lon)
	bathy = np.array(bathy)

	d = np.tile(bathy.min()-10,bathy.shape)


	X,Y = np.meshgrid(adcp.columns.values,-adcp.index.values)
	vlevs = np.arange(-1.5,1.5,0.1)

	fig,ax = plt.subplots(figsize=(15,7))
	ax.contourf(X,Y,adcp,100,cmap='RdBu_r',vmin=-1.5,vmax=1.5)
	c1 = ax.contour(X,Y,adcp,levels=vlevs,colors='k')
	ax.clabel(c1,fmt='%1.1f', fontsize=10)

	X,Y = np.meshgrid(dens0.columns.values,-dens0.index.values)
	dlevs = np.array([24.5,26.75,27])
	
	c2 = ax.contour(X,Y,dens0,levels=dlevs,colors='k',linewidths=2,linestyles='dashed')
	clabels = ax.clabel(c2, fontsize=10,
	fmt='%1.2f',colors='k',use_clabeltext=True)
	
	plt.ylim([lim,0])
	

	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))
		
	ax.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)
	
	
	ax.plot(adcp.columns.values,np.repeat(0,adcp.columns.values.shape),
	marker='v',markersize=20,color='k',
	linewidth=0)



def plotdf(secvel,df,lim=-400):
	adcp = df
	dens0 = secvel.dens0
	
	lon = secvel.lon.mean().values
	lat = secvel.lat.mean().values

	a,bathy = etopo_sec(lat,lon)
	bathy = np.array(bathy)

	d = np.tile(bathy.min()-10,bathy.shape)


	X,Y = np.meshgrid(adcp.columns.values,-adcp.index.values)
	vlevs = np.arange(-17,17,1)

	fig,ax = plt.subplots(figsize=(15,7))
	ax.contourf(X,Y,adcp,100,cmap='RdBu_r',vmin=-17,vmax=17)
	c1 = ax.contour(X,Y,adcp,levels=vlevs,colors='k')
	ax.clabel(c1,fmt='%1.1f', fontsize=10)

	X,Y = np.meshgrid(dens0.columns.values,-dens0.index.values)
	dlevs = np.array([24.5,26.5,26.75,27])
	
	c2 = ax.contour(X,Y,dens0,levels=dlevs,colors='k',linewidths=2,linestyles='dashed')
	clabels = ax.clabel(c2, fontsize=10,fmt='%1.1f',colors='k',use_clabeltext=True)
	
	plt.ylim([lim,0])
	

	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))
		
	ax.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)
	
	ax.plot(adcp.columns.values,np.repeat(0,adcp.columns.values.shape),
	marker='v',markersize=20,color='k',
	linewidth=0)
	
	
	

def plotsec(opname,rad,pastafig,
	pathadcp,corsta='b',cond=False,
	ts=False,fmed=False,interp=False,
	radcomplicada=False,mdrcalc=True):


	vsec,latv,lonv,lat,lon,section,lonsta,latsta = [],[],[],[],[],[],[],[]
	adcp,lat_adcp,lon_adcp = [],[],[]
	
	radname = 'r'+os.path.split(rad)[1]
	
	if mdrcalc:
		vsec,lat,lon,section = extsect(rad,sup=True)
	else:	
		vsec,lat,lon,section = extsect(rad,sup=False)
	

	
#	opname = section.items.values[0].split('0')[0]

	if ts:
		tscomp(section)
	
		plt.savefig(pastafig+opname+'/ts/'+'ts'+radname+'.png',dpi=100,format='png')
		plt.close('all')	

	lat_adcp,lon_adcp,adcp = adcpvel(lat,lon,pathadcp,dya=600,interp=interp)

	if mdrcalc:
		vsec = mdr(vsec,adcp)
	
#	lonsta = np.nanmean(section.minor_xs('lon'),axis=0)
#	latsta = np.nanmean(section.minor_xs('lat'),axis=0)
	
	latv = lat[:-1]+np.diff(lat)/2
	lonv = lon[:-1]+np.diff(lon)/2
	
	if fmed:
		f = sw.f(latv)
		vsec = vsec*(f/np.mean(f))
	
	
	if cond:
	#	SALINIDADE
		plotprop(z = section.minor_xs('sp'),
		lat=latsta,lon=lonsta,cor=corsta,csteps=0.5,
		title=u'Salinidade Prática',
		propmin=34,propmax=37)

		plt.savefig(pastafig+opname+'/salinidade/'+'sal'+radname+'.png',dpi=100,format='png')
		plt.close('all')

	#	TEMPERATURA	
		plotprop(z = section.minor_xs('pt'),
		lat=latsta,lon=lonsta,cor=corsta,csteps=2.5,
		title=u'Temperatura Potencial ($\\theta_0$)',
		propmin=0,propmax=28)
	
		plt.savefig(pastafig+opname+'/temperatura/'+'temp'+radname+'.png',dpi=100,format='png')
		plt.close('all')

	#	MASSAS D'ÁGUA
		plotprop(z = section.minor_xs('psigma0'),
		lat=latsta,lon=lonsta,cor=corsta,wmass=True,
		title=u'Massas D\'água ($\\rho_0$)')
	
		plt.savefig(pastafig+opname+'/massas/'+'mass'+radname+'.png',dpi=100,format='png')
		plt.close('all')
	
	#	DENSIDADE POTENCIAL	
		plotprop(z = section.minor_xs('psigma0'),
		lat=latsta,lon=lonsta,cor=corsta,csteps=0.5,
		title=u'Densidade Potencial ($\\rho_0$)',
		propmin=23,propmax=28)
		
		plt.savefig(pastafig+opname+'/densidade/'+'dens'+radname+'.png',dpi=100,format='png')
		plt.close('all')

	
#	reff = adcp.iloc[150,:].values-vsec.iloc[150,:].values
#	vsec += reff
	if radcomplicada:
		latv[3] = -2.0933934592636545
		lonv[3] = -38.442823159524877

	if mdrcalc:
		ttitle = u'Velocidade Geostrófica (MDR)'
	else:
		ttitle = u'Velocidade Geostrófica (Isopicnal = 32.15)'
		
#	Velocidade Geostrófica	
	plotprop(z = vsec,cmap = 'RdBu_r',
	lat=latv,lon=lonv,cor=corsta,csteps=0.1,
	title=ttitle,
	propmin=-1.5,propmax=1.5,lim=-300)
		
	plt.savefig(pastafig+opname+'/geostrofia/'+'geos'+radname+'.png',dpi=100,format='png')
	plt.close('all')

	plotprop(z = vsec,lat=latv,lon=lonv,cor=corsta,
	z2 = adcp,lat2=lat_adcp,lon2=lon_adcp,cor2='g',
	csteps=0.1,title=u'Velocidade Geostrófica X Velocidade Observada',
	cmap = 'RdBu_r',propmin=-1.5,propmax=1.5,comp=True,lim=-300)
	
	plt.savefig(pastafig+opname+'/velocidade/'+'vel'+radname+'adcp.png',dpi=100,format='png')
	plt.close('all')


def plotvel(secvel,velname='velmdr',
	lim=-400,radcomplicada=False,
	figsize = (15,7),llcrnrlon=-51,
	urcrnrlon=-30.2, llcrnrlat=-15,
	urcrnrlat=7.1,projection='merc',
	resolution='l',cor='r',refadcp=None,
	mmax=1.3):


	cmapp = 'RdBu_r'
	alp = 1
#	mmax = 1.3


	fig = plt.figure(figsize=figsize)

	exec('adcp = secvel.'+velname)
	dens0 = secvel.dens0
	dens1 = secvel.dens1
	dens2 = secvel.dens2
	dens4 = secvel.dens4

	
	nnan = adcp.columns.values[~np.isnan(adcp.mean().values)]
	
	lon = secvel.lon.mean().values
	lat = secvel.lat.mean().values

	if radcomplicada:
		lat[3] = -2.0933934592636545
		lon[3] = -38.442823159524877

	a,bathy = etopo_sec(lat,lon)
	bathy = np.array(bathy)

	d = np.tile(bathy.min()-10,bathy.shape)


	X,Y = np.meshgrid(adcp.columns.values,-adcp.index.values)
	vlevs = np.arange(-mmax,mmax+0.1,0.1)

	gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1])

#	fig,ax = plt.subplots(figsize=(15,7))
	
	ax1 = plt.subplot(gs[0,:])
	
#	CF = ax1.contourf(X,Y,adcp,vlevs,cmap='RdBu_r',vmin=-1.5,vmax=1.5)
	CF = ax1.contourf(X,Y,adcp,vlevs,cmap=cmapp,alpha= alp,
		vmin=-mmax,vmax=mmax)
	c1 = ax1.contour(X,Y,adcp,levels=vlevs,colors='k')
	ax1.clabel(c1,fmt='%1.1f', fontsize=10)

	try:
		ax1.plot(adcp.columns.values.astype('float'),
		refadcp,color='r',linewidth=3)
	except:
		print 'Sem nível de referência.'

	X,Y = np.meshgrid(dens0.columns.values,-dens0.index.values)


	dlevs = np.array([24.5,26.8])
	
	c2 = ax1.contour(X,Y,dens0,levels=dlevs,colors='k',
	linewidths=2,linestyles='dashed')
	clabels = ax1.clabel(c2, fontsize=10,
	fmt='%1.2f',colors='k',use_clabeltext=True)
	
	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))	



#	dlevs = np.array([32.15])
#	
#	c2 = ax1.contour(X,Y,dens1,levels=dlevs,colors='k',
#	linewidths=2,linestyles='dashed')
#	clabels = ax1.clabel(c2, fontsize=10,
#	fmt='%1.2f',colors='k',use_clabeltext=True)
#	
#	for txt in clabels:
#		txt.set_rotation(0)
#		txt.set_bbox(dict(facecolor='white',
#		alpha=0.8,edgecolor='None'))

#	dlevs = np.array([37])

#	c2 = ax1.contour(X,Y,dens2,levels=dlevs,colors='k',
#	linewidths=2,linestyles='dashed')
#	clabels = ax1.clabel(c2, fontsize=10,
#	fmt='%1.2f',colors='k',use_clabeltext=True)
#	
#	for txt in clabels:
#		txt.set_rotation(0)
#		txt.set_bbox(dict(facecolor='white',
#		alpha=0.8,edgecolor='None'))
#			
#	dlevs = np.array([45.83,45.80])

#	c2 = ax1.contour(X,Y,dens4,levels=dlevs,colors='k',
#	linewidths=2,linestyles='dashed')
#	clabels = ax1.clabel(c2, fontsize=10,
#	fmt='%1.2f',colors='k',use_clabeltext=True)
#	
#	for txt in clabels:
#		txt.set_rotation(0)
#		txt.set_bbox(dict(facecolor='white',
#		alpha=0.8,edgecolor='None'))

	ax1.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)

	ax1.set_ylim([lim,0])

	ax1.set_xlim([nnan[0],nnan[-1]])



	
	ax1.plot(adcp.columns.values,np.repeat(0,adcp.columns.values.shape),
	marker='v',markersize=20,color='k',
	linewidth=0)


	ax2 = plt.subplot(gs[1,:])
	
#	ax2.contourf(X,Y,adcp,vlevs,cmap='RdBu_r',vmin=-1.5,vmax=1.5)
	ax2.contourf(X,Y,adcp,vlevs,cmap=cmapp,alpha=alp,
		vmin=-mmax,vmax=mmax)
	c1 = ax2.contour(X,Y,adcp,levels=vlevs,colors='k')
	ax2.clabel(c1,fmt='%1.1f', fontsize=10)

	X,Y = np.meshgrid(dens0.columns.values,-dens0.index.values)


#	dlevs = np.array([24.5,26.8])
#	
#	c2 = ax2.contour(X,Y,dens0,levels=dlevs,colors='k',
#	linewidths=2,linestyles='dashed')
#	clabels = ax2.clabel(c2, fontsize=10,
#	fmt='%1.2f',colors='k',use_clabeltext=True)
#	
#	for txt in clabels:
#		txt.set_rotation(0)
#		txt.set_bbox(dict(facecolor='white',
#		alpha=0.8,edgecolor='None'))	

	dlevs = np.array([32.15])
	
	c2 = ax2.contour(X,Y,dens1,levels=dlevs,colors='k',
	linewidths=2,linestyles='dashed')
	clabels = ax2.clabel(c2, fontsize=10,
	fmt='%1.2f',colors='k',use_clabeltext=True)
	
	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))

	dlevs = np.array([37])

	c2 = ax2.contour(X,Y,dens2,levels=dlevs,colors='k',
	linewidths=2,linestyles='dashed')
	clabels = ax2.clabel(c2, fontsize=10,
	fmt='%1.2f',colors='k',use_clabeltext=True)
	
	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))
			
	dlevs = np.array([45.83,45.80])

	c2 = ax2.contour(X,Y,dens4,levels=dlevs,colors='k',
	linewidths=2,linestyles='dashed')
	clabels = ax2.clabel(c2, fontsize=10,
	fmt='%1.2f',colors='k',use_clabeltext=True)
	
	for txt in clabels:
		txt.set_rotation(0)
		txt.set_bbox(dict(facecolor='white',
		alpha=0.8,edgecolor='None'))


	ax2.set_ylim([np.min(np.min(Y)),lim])
	
	ax2.set_xlim([nnan[0],nnan[-1]])

	ax2.fill_between(adcp.columns.values.astype('float'),
	bathy, d, interpolate=True, color='black',zorder=20)
	
	ax2.plot(adcp.columns.values,np.repeat(0,adcp.columns.values.shape),
	marker='v',markersize=20,color='k',
	linewidth=0)
	
	position=fig.add_axes([0.91,0.165,0.01,0.735])
#	fig.colorbar(CF,orientation='horizontal',
#	cax=position,extend='both')
	
	cbar = fig.colorbar(CF,
		ticks=np.arange(-mmax,mmax+0.1,0.2),cax=position,
		extend='both')
	
	pos = np.array(ax1.get_position().bounds)
	pos[1] -= pos[1]*0.69
	ax2.set_position(pos)
	plt.draw()
	

#	### MAP PLOT###
#	# starts with the second subplot
#	plt.subplot(gs[:,1])
#	# create a map with m as handle with attributes
#	# defined in the begin

#	m = Basemap(llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
#                llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
#                projection=projection, resolution=resolution)


#	# draw states and coastlines as default
#	m.drawstates()
#	m.drawcoastlines()
#	# fill continents with gray color
#	m.fillcontinents(color='0.85')
#	# create meridians array that will be from ceil of the low limit
#	# to floor of the up limit plus 2 with 5 steps
#	meridians = np.arange(np.ceil(llcrnrlon), np.floor(urcrnrlon) + 2, 5)
#	# create the parallels array that will be the round of the 10 values
#	# between the floor of the low limit to floor of up limit
#	parallels = np.round(np.linspace(np.floor(llcrnrlat),np.floor(urcrnrlat),10))
#	# draw parallels and meridians with 0.4 of line and this pattern of labels
#	m.drawparallels(parallels, linewidth=0.4, labels=[1, 0, 0, 0]) # labels on the left
#	m.drawmeridians(meridians, linewidth=0.4, labels=[0, 0, 0, 1]) # labels on the bottom

#	plt.title(u'Localização')	

#	lon = secvel.lon.mean().values
#	lat = secvel.lat.mean().values

#	# find the coordiates for the median values of latitudes and longitudes
#	lonii,latii = m(lon,lat)
#	# plot it on the map with defined color as circles with no line between their
#	m.plot(lonii,latii,color=cor,alpha=0.9,markersize=4, 
#		marker='o',antialiased=True, linewidth=0)
#		
#				


#	# make etopo background with 0.3 of alpha
#	m.drawmapboundary(fill_color='#D1F0FF')


	fig.text(0.08,0.5,'Profundidade',rotation=90,
		ha='center',va='center')
	ax2.set_xlabel(u'Distância',fontweight='bold')

#	return fig,ax1,ax2,m
	return fig,ax1,ax2



	
	
	
	
#def plotprop(z,lat,lon,lim=-600, figsize = (26,6),
#	llcrnrlon=-51, urcrnrlon=-30.2, llcrnrlat=-15,
#	urcrnrlat=7.1,projection='merc',resolution='l',
#	cor='r',csteps=0.5,cmap='Spectral_r',
#	title='Titulo',wmass=False,etopoflag=False,
#	propmin=[],propmax=[],comp=False,z2=[],
#	lat2=[],lon2=[],cor2=[],erro=[],err=False):
#	
#	'''
#	This function plot some property section
#	with its path in map
#	'''
#	
#	
#	
#	# se temos nan à mil metros da primeira estação
#	# usa-se isnull pq o nan do DataFrame não é reconhecido
#	# pelo numpy
##	if pd.isnull(z.values[1800,-1]):
##		print ' A Radial está invertida! '
##		print ' --- --- invertendo a radial'
##		z = z.reindex(columns=np.flipud(z.columns))
##		lat = np.flipud(lat)
##		lon = np.flipud(lon)
##		print '###########################'
#	
#	
#	dist,bathy = etopo_sec(lat,lon)
#	dist = np.array(dist)
#	bathy = np.array(bathy)
#	
#	x = dist
#	y = -np.array(z.index.values)
#	
#	if comp:
#		x2 = np.hstack([0,np.cumsum(sw.dist(lat2,lon2)[0])])
#		y2 = -np.array(z2.index.values)
#	

#	### PREPARING THE FIGURE ###
#	# create a figure with figsize
#	fig = plt.figure(figsize=figsize)
#	# create a gridspec for width ratios in subplots
#	# this is used because we want that the figure have 2
#	# subplots in line and the first must be the double of
#	# the second's width
#	if err==True:
#		gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1])
#	else:
#		gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1])
#	

#	if not propmin:
#		propmin = np.min(np.min(z))
#		propmax = np.max(np.max(z))
#	
#	# select the first axis
#	ax1 = plt.subplot(gs[0,0])


#	if wmass:
#		### CONTOUR PLOT ###
#		# determine the velocity contour levels
#		clev = [20,25.6,26.9,27.38,27.53,28]
#		norml = BoundaryNorm(clev, 256)
#		
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		ax1.contourf(x,y,z,levels=clev,alpha=0.5,
#		colors=('#ff0000', '#ff9900', '#00c4ff', '#4fda66', '#0000ff'),
#		norm = norml)
#	else:
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		ax1.contourf(x,y,z,100,
#		cmap=cmap,vmin=propmin,vmax=propmax)

#		### CONTOUR PLOT ###
#		# determine the velocity contour levels
#		clev = np.arange(propmin,propmax,csteps)

#	# make contour plot with vlev levels in black
#	c2 = ax1.contour(x,y,z,levels = clev, colors='k')


#	# make contour lables with 10 size
#	plt.clabel(c2, fontsize=10)

######	batimetria	

#	d = np.tile(bathy.min()-10,bathy.shape)
#	ax1.fill_between(dist, bathy, d,
#	interpolate=True, color='black',alpha=0.6)
#			
#	### CHANGING THE AXIS LIMITS ###
#	# set the y lim from lim variable (defalult is -600)
#	ax1.set_ylim([lim,0])
#	# set xlim to hide nan columns from plot (read extvel to know
#	# why this happens)
########	plt.xlim([x[~np.isnan(z[10,:])][0],
########		  x[~np.isnan(z[10,:])][-1]])
#	# create a grid on plot from the ticks
#	ax1.grid('on')
#	plt.setp(ax1.get_xticklabels(), visible=False)

#	if comp:
#		ax1.plot(x,np.repeat(lim,x.shape),
#		marker='^',markersize=20,color='k',
#		linewidth=0)
#		plt.draw()

#	ax2 = plt.subplot(gs[1,0])
#	
#	if wmass:
#		### CONTOUR PLOT ###
#		# determine the velocity contour levels
#		clev = [20,25.6,26.9,27.38,27.53,28]
#		norml = BoundaryNorm(clev, 256)
#		
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		c1 = ax2.contourf(x,y,z,levels=clev,alpha=0.5,
#		colors=('#ff0000', '#ff9900', '#00c4ff', '#4fda66', '#0000ff'),
#		norm = norml)
#		c2 = ax2.contour(x,y,z,levels = clev, colors='k')
#	elif comp:
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		c1 = ax2.contourf(x2,y2,z2,100,
#		cmap=cmap,vmin=propmin,vmax=propmax)

#		### CONTOUR PLOT ###
#		# determine the velocity contour levels
#		clev = np.arange(propmin,propmax,csteps)
#		# make contour plot with vlev levels in black
#		c2 = ax2.contour(x2,y2,z2,levels = clev, colors='k')
#		
#	else:
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		c1 = ax2.contourf(x,y,z,100,
#		cmap=cmap,vmin=propmin,vmax=propmax)

#		### CONTOUR PLOT ###
#		# determine the velocity contour levels
#		clev = np.arange(propmin,propmax,csteps)
#	
#		# make contour plot with vlev levels in black
#		c2 = ax2.contour(x,y,z,levels = clev, colors='k')

#	# make contour lables with 10 size
#	plt.clabel(c2, fontsize=10)
#	
#	# plot the colorbar with cbar as handle
#	position=fig.add_axes([0.12,0.06,0.48,0.02])
#	fig.colorbar(c1,orientation='horizontal',cax=position)

#	ax2.fill_between(dist, bathy, d,
#	interpolate=True, color='black',alpha=0.6)

#	if comp:
#	
#		### CHANGING THE AXIS LIMITS ###
#		# set the y lim from lim variable (defalult is -600)
#		ax2.set_ylim([lim,0])
#		# set xlim to hide nan columns from plot (read extvel to know
#		# why this happens)
#		# create black bottom triangles to show where we have profile data
#		ax2.plot(x2,np.repeat(lim,x2.shape),
#		marker='^',markersize=20,color='k',
#		linewidth=0)
#		plt.draw()
#			
#	else:
#	
#		### CHANGING THE AXIS LIMITS ###
#		# set the y lim from lim variable (defalult is -600)
#		ax2.set_ylim([np.min(np.min(y)),lim])
#		# set xlim to hide nan columns from plot (read extvel to know
#		# why this happens)
#		# create black bottom triangles to show where we have profile data
#		ax2.plot(x,np.repeat(np.min(np.min(y)),x.shape),
#		marker='^',markersize=20,color='k',
#		linewidth=0)

#	#######	plt.xlim([x[~np.isnan(z[10,:])][0],
#	#######		  x[~np.isnan(z[10,:])][-1]])

#	# create a grid on plot from the ticks
#	ax2.grid('on')

#	
#	pos = np.array(ax1.get_position().bounds)
#	pos[1] -= pos[1]*0.69
#	ax2.set_position(pos)
#	plt.draw()
#	
#		
#	### SETTING LABELS ###
#	# make contour lables with 10 size
#	plt.clabel(c2, fontsize=10)
#	# make x and y labels
#	plt.figtext(figure=fig,x=0.09,y=0.5,s=u'Profundidade (m)',rotation=90,va='center')
#	ax2.set_xlabel(u'Distância (Km)')

#	# make title, $ is used to wite in latex and u'' is used 
#	# to show that its wrote in utf-8
#	ax1.set_title(title)

#	if wmass:
#		position.text(x=0.1,y=-1.5,ha='center',s='AT')
#		position.text(x=0.3,y=-1.5,ha='center',s='ACAS')
#		position.text(x=0.5,y=-1.5,ha='center',s='AIA')
#		position.text(x=0.7,y=-1.5,ha='center',s='ACS')
#		position.text(x=0.9,y=-1.5,ha='center',s='APAN')

#	if err:
#		ax3 = plt.subplot(gs[2,0])
#		### FILLED CONTOUR ###
#		# plot the filled contour of the data with 100 colors
#		# inverted RedtoBlue cmap and colorbar limit from -2 to 2
#		c3 = ax3.contourf(x2,y2,erro,100,
#		vmin=0,vmax=1)

##		### CONTOUR PLOT ###
##		# determine the velocity contour levels
##		clev3 = np.arange(0,1,0.2)
##		# make contour plot with vlev levels in black
##		c3_2 = ax3.contour(x2,y2,erro,levels = clev3, colors='k')
##		plt.clabel(c3, fontsize=10)

#	### MAP PLOT###
#	# starts with the second subplot
#	plt.subplot(gs[:,1])
#	# create a map with m as handle with attributes
#	# defined in the begin

#	m = Basemap(llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
#                llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
#                projection=projection, resolution=resolution)


#	# draw states and coastlines as default
#	m.drawstates()
#	m.drawcoastlines()
#	# fill continents with gray color
#	m.fillcontinents(color='0.85')
#	# create meridians array that will be from ceil of the low limit
#	# to floor of the up limit plus 2 with 5 steps
#	meridians = np.arange(np.ceil(llcrnrlon), np.floor(urcrnrlon) + 2, 5)
#	# create the parallels array that will be the round of the 10 values
#	# between the floor of the low limit to floor of up limit
#	parallels = np.round(np.linspace(np.floor(llcrnrlat),np.floor(urcrnrlat),10))
#	# draw parallels and meridians with 0.4 of line and this pattern of labels
#	m.drawparallels(parallels, linewidth=0.4, labels=[1, 0, 0, 0]) # labels on the left
#	m.drawmeridians(meridians, linewidth=0.4, labels=[0, 0, 0, 1]) # labels on the bottom

#	plt.title(u'Localização')	

#	if comp:
#		# find the coordiates for the median values of latitudes and longitudes
#		lonii,latii = m(lon2,lat2)
#		# plot it on the map with defined color as circles with no line between their
#		m.plot(lonii,latii,color=cor2,alpha=0.9,markersize=4, 
#			marker='o',antialiased=True, linewidth=0)
#		
#				
#	# find the coordiates for the median values of latitudes and longitudes
#	loni,lati = m(lon,lat)
#	# plot it on the map with defined color as circles with no line between their
#	m.plot(loni,lati,color=cor,alpha=0.9,markersize=6, 
#		marker='o',antialiased=True, linewidth=0)

#	# make etopo background with 0.3 of alpha
#	if etopoflag:
#		m.etopo(alpha=0.3)
#	else:
#		m.drawmapboundary(fill_color='#D1F0FF')


#	return fig,ax1,ax2,m

