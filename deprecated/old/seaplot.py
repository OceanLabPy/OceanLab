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

#from pandas import Panel

import numpy as np


#import scipy.interpolate as interp
import scipy.io as sio
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

#import fnmatch


import gsw
import seawater as sw
import matplotlib.patches
import matplotlib


	
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
    

def tscomp(secref,path='/home/iury/Copy/TCC/dados/CLIMAT/WOA'):
	'''
	This function creates T-S plots for CTD stations
	and referred climatology from WOA (NOAA).
	'''
	

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
	
