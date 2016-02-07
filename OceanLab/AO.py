# -*- coding: utf-8 -*-
from scipy.optimize import curve_fit
import itertools
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def psi2uv(x,y,psi):
    '''
    % PSI2UV Velocity components from streamfunction
    %   [U,V] = PSI2UV(X,Y,PSI) returns the velocity components U and V
    %   from streamfunction PSI. X and Y are longitude and latitude matrices,
    %   respectively, as PSI. All matrices have same dimension.
    % 
    %   --> !!! IMPORTANT !!!  <-----------------------------------------
    %   The grid indexes IM and JM must have origin on the left-inferior 
    %   corner, increasing to right and to up, respectively.
    % 
    %                              GRID
    %            :            :            :            :
    %            |            |            |            |
    %   (JM = 2) o ---------- o ---------- o ---------- o --- ...
    %            |            |            |            |
    %            |            |            |            |
    %   (JM = 1) o ---------- o ---------- o ---------- o --- ...
    %            |            |            |            |
    %            |            |            |            |
    %   (JM = 0) o ---------- o ---------- o ---------- o --- ...
    %        (IM = 0)     (IM = 1)     (IM = 2)     (IM = 3)
    % 
    %   -----------------------------------------------------------------
    % 
    %   USAGE:  u,v = psi2uv(x,y,psi)
    % 
    %   INPUT:
    %     x     = decimal degrees (+ve E, -ve W) [-180..+180]
    %     y     = decimal degrees (+ve N, -ve S) [- 90.. +90]
    %     psi   = streamfunction [m^2 s^-1]
    % 
    %   OUTPUT:
    %     u   = velocity zonal component [m s^-1]
    %     v   = velocity meridional component [m s^-1]
    % 
    %   EXAMPLE:
      import numpy as np
      import matplotlib.pyplot as plt
      x,y = np.meshgrid(np.arange(-30,-20,0.5),np.arange(-35,-25,0.5))
      psi = (x-np.mean(np.mean(x)))**2 + (y-np.mean(np.mean(y)))**2
    
      plt.ion()
      plt.figure()
      u,v = psi2uv(x,y,psi)
      plt.contourf(x,y,psi,30,cmap='Spectral_r')
      plt.quiver(x,y,u,v,color='w')

    %   Author:
    %   Rafael A. de Mattos (RAM), 03 May 2005
    %   Update of 03 May 2005 (RAM)
       Translated to PYTHON by Iury Sousa, 29 July 2015

    %   ======================================================================
    '''

    JM,IM = psi.shape

    u = np.zeros((JM,IM))
    v = u
    
    u2,v2 = u.copy()*np.nan,v.copy()*np.nan
    u1,v1 = u.copy()*np.nan,v.copy()*np.nan

    #% Velocity components u2 and v2 on Natural Coordinates

    for i in np.arange(0,IM):
     for j in np.arange(0,JM):
    
        #% Take for behind difference
        if j == 0:
            dy2 = (y[j+1,i]-y[j,i])*60*1852
            dx2 = (x[j+1,i]-x[j,i])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j,i]-psi[j+1,i])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90

        v2[j,i] = V*np.sin(ang*np.pi/180)
        u2[j,i] = V*np.cos(ang*np.pi/180)

        
        
        #% Take forward differences
        if j == JM-1:
            dy2 = (y[j,i]-y[j-1,i])*60*1852
            dx2 = (x[j,i]-x[j-1,i])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j-1,i]-psi[j,i])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90

        v2[j,i] = V*np.sin(ang*np.pi/180)
        u2[j,i] = V*np.cos(ang*np.pi/180)

        
        #% Take centered differences on interior points
        if (j > 0) & (j < JM-1):
            dy2 = (y[j+1,i]-y[j-1,i])*60*1852
            dx2 = (x[j+1,i]-x[j-1,i])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j-1,i]-psi[j+1,i])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90
        v2[j,i] = V*np.sin(ang*np.pi/180)
        u2[j,i] = V*np.cos(ang*np.pi/180)


    #% Velocity components u1 and v1 on Natural Coordinates

    for j in np.arange(0,JM):
     for i in np.arange(0,IM):

        #% Take for behind difference
        if i == 0:
            dy2 = (y[j,i+1]-y[j,i])*60*1852
            dx2 = (x[j,i+1]-x[j,i])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j,i]-psi[j,i+1])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90

        v1[j,i] = V*np.sin(ang*np.pi/180)
        u1[j,i] = V*np.cos(ang*np.pi/180)

        
        
        #% Take forward differences
        if i == IM-1:
            dy2 = (y[j,i]-y[j,i-1])*60*1852
            dx2 = (x[j,i]-x[j,i-1])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j,i-1]-psi[j,i])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90

        v1[j,i] = V*np.sin(ang*np.pi/180)
        u1[j,i] = V*np.cos(ang*np.pi/180)

        
        #% Take centered differences on interior points
        if (i > 0) & (i < IM-1):
            dy2 = (y[j,i+1]-y[j,i-1])*60*1852
            dx2 = (x[j,i+1]-x[j,i-1])*60*1852
            dd = np.sqrt(dx2**2+dy2**2)
            V = (psi[j,i-1]-psi[j,i+1])/dd
        if np.isinf(dy2/dx2):
            ang = 0
        elif (np.math.atan(dy2/dx2)*180/np.pi > 0) & (~np.isinf(dy2/dx2)):
            ang = np.math.atan(dy2/dx2)*180/np.pi - 90
        else:
            ang = np.math.atan(dy2/dx2)*180/np.pi + 90
        v1[j,i] = V*np.sin(ang*np.pi/180)
        u1[j,i] = V*np.cos(ang*np.pi/180)
    
    #% Components u and v on Cartesian Coordinates

    v = v2-v1
    u = u2-u1
    
    #on original the final was
    #v = v1+v2
    #u = u1+u2

    return u,v

def vectoa(xc,yc,x,y,u,v,corrlenx,corrleny,err,b):
        '''
	#%function psi=vectoa(xc,yc,x,y,u,v,corrlenx,corrleny,err,b)
	#%       (xc,yc) are vectors (row or column) of interpolation points
	#%       (x,y)   are vectors (row or column) of observation points
	#%       (u,v)   are matrices of east and north components with each day
	#%               entered column-wise or row-wise
	#%       corrlen,err are correlation length scales (in x and y) and error for a
	#%               gaussian streamfunction covariance function
	#%       psi     is streamfunction at    (xc,yc) returned in same
	#%               format as u or v
	#PYTHON VERSION MADE BY Iury Simões - Oceanography from Universidade Federal do Ceará
	# 2014
        '''
        corrlen = corrleny
        xc *= corrleny*1./corrlenx
        x *= corrleny*1./corrlenx
        
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
		t.append(np.math.atan2(ii,jj))
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
		tc.append(np.math.atan2(ii,jj))
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


######################################################################
##########################SCALOA######################################
######################################################################
def scaloa(xc, yc, x, y, t=None, corrlenx=None,corrleny=None, err=None, zc=None):
  """
  Scalar objective analysis. Interpolates t(x, y) into tp(xc, yc)
  Assumes spatial correlation function to be isotropic and Gaussian in the
  form of: C = (1 - err) * np.exp(-d**2 / corrlen**2) where:
  d : Radial distance from the observations.
  Parameters
  ----------
  corrlen : float
  Correlation length.
  err : float
  Random error variance (epsilon in the papers).
  Return
  ------
  tp : array
  Gridded observations.
  ep : array
  Normalized mean error.
  Examples
  --------
  See https://ocefpaf.github.io/python4oceanographers/blog/2014/10/27/OI/
  Notes
  -----
  The funcion `scaloa` assumes that the user knows `err` and `corrlen` or
  that these parameters where chosen arbitrary. The usual guess are the
  first baroclinic Rossby radius for `corrlen` and 0.1 e 0.2 to the sampling
  error.
  """
  corrlen = corrleny
  xc *= corrleny*1./corrlenx
  x *= corrleny*1./corrlenx
  
  n = len(x)
  x, y = np.reshape(x, (1, n)), np.reshape(y, (1, n))
  # Squared distance matrix between the observations.
  d2 = ((np.tile(x, (n, 1)).T - np.tile(x, (n, 1))) ** 2 +
  (np.tile(y, (n, 1)).T - np.tile(y, (n, 1))) ** 2)
  nv = len(xc)
  xc, yc = np.reshape(xc, (1, nv)), np.reshape(yc, (1, nv))
  # Squared distance between the observations and the grid points.
  dc2 = ((np.tile(xc, (n, 1)).T - np.tile(x, (nv, 1))) ** 2 +
  (np.tile(yc, (n, 1)).T - np.tile(y, (nv, 1))) ** 2)
  # Correlation matrix between stations (A) and cross correlation (stations
  # and grid points (C))
  A = (1 - err) * np.exp(-d2 / corrlen ** 2)
  C = (1 - err) * np.exp(-dc2 / corrlen ** 2)
  if 0: # NOTE: If the parameter zc is used (`scaloa2.m`)
    A = (1 - d2 / zc ** 2) * np.exp(-d2 / corrlen ** 2)
    C = (1 - dc2 / zc ** 2) * np.exp(-dc2 / corrlen ** 2)
  # Add the diagonal matrix associated with the sampling error. We use the
  # diagonal because the error is assumed to be random. This means it just
  # correlates with itself at the same place.
  A = A + err * np.eye(len(A))
  # Gauss-Markov to get the weights that minimize the variance (OI).
  tp = None
  ep = 1 - np.sum(C.T * np.linalg.solve(A, C.T), axis=0) / (1 - err)
  if t!=None:
    t = np.reshape(t, (n, 1))
    tp = np.dot(C, np.linalg.solve(A, t))
    #if 0: # NOTE: `scaloa2.m`
    #  mD = (np.sum(np.linalg.solve(A, t)) /
    #  np.sum(np.sum(np.linalg.inv(A))))
    #  t = t - mD
    #  tp = (C * (np.linalg.solve(A, t)))
    #  tp = tp + mD * np.ones(tp.shape)
    return tp, ep
    
  if t==None:
    print("Computing just the interpolation errors.")
    #Normalized mean error. Taking the squared root you can get the
    #interpolation error in percentage.
    return ep
######################################################################
##########################SCALOA_end######################################
######################################################################



def combinations(array):
	'''
	Calcula todas as combinacoes possiveis
	de uma array com ela mesma
	'''
	comb = []
	for n in itertools.combinations(array,2):
	  comb.append(n)
	return np.array(comb)

def cor_calc(cdat,dist,binn,cut=20):
	dist = np.squeeze(dist)
	lims = np.arange(0,np.ceil(dist.max())+binn,binn)
	rs = lims[:-1]+np.diff(lims)/2
	
	cor = []
	bin_pt = []
	for inf,sup in zip(lims[:-1],lims[1:]):
		args = np.argwhere((dist>inf)&(dist<=sup))
		bin_pt.append(args.size)
		if (args.size==0)|(args.size==1):
		#if (args.size<cut):
			cor.append(np.nan)
		else:
		        cc = np.squeeze(cdat[args,:])
		        cc = (cc-cc.mean(axis=0))
		        CP = cc[:,0]*cc[:,1]
		        cor.append(np.mean(CP)/(cc.std(axis=0).prod()))
		        
	cor = np.array(cor)


	return cor,rs,bin_pt

#FUNÇÕES PARA O AJUSTE

def func(x,A,lc):
	'''
	funcao de ajuste da correlacao
	'''
	return A*np.exp((-(x)**2)/(lc**2))
    
def adjusting(cdat,dist,binn,cut=550):
    corr,rrs,b_pt = cor_calc(cdat,dist,binn)
    
    ydata = corr[(~np.isnan(corr))&(rrs<cut)]
    xdata = rrs[(~np.isnan(corr))&(rrs<cut)]

    xdata = np.hstack([-np.flipud(xdata),xdata])
    ydata = np.hstack([np.flipud(ydata),ydata])

    #f = scint.interp1d(xdata,ydata)
    #X = np.arange(xdata.min(),xdata.max())
    #Y = f(X)
    
    X = xdata
    Y = ydata

    popt1,popv = curve_fit(func,X,Y,p0=(1,binn),maxfev=100000)

    LIM = cut
    yy = func(np.arange(-LIM,LIM),*popt1)

    #plt.close()
    plt.figure()
    plt.scatter(xdata,ydata)
    plt.plot(np.arange(-LIM,LIM),yy,'r',linewidth=2)
    plt.ylim([-0.1,1.2])
    plt.xlim([0,LIM])
    plt.grid('on')
    plt.xlabel('Lag [Km]')
    plt.ylabel(r'Correla\c{c}\~ao')
    
    ER = 1-func(0,*popt1)
    LC = popt1[-1]
    plt.text(LIM*0.5,0.4,
            r'$\epsilon^{2}$: %.3f - LC: %.3f km'%(ER,LC),
            fontsize=15)
    
    return ER,LC


#BOUNDARY CONDITION

def bathy_lims(PROF,lnu=-33,lnd=-49,
                ltu=-20,ltd=-34,step=1,
                etopo=5,just_bathy=False):
    
    arq = 'etopo%.i.nc'%(etopo)
    #LENDO DADOS DO ETOPO
    url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/'+arq
    etopodata = Dataset(url)
    
    if etopo==1:
        strx = 'x'
        stry = 'y'
        strb = 'rose'
        blon = (etopodata.variables[strx][:])
        blat = etopodata.variables[stry][:]
    elif etopo==5:
        strx = 'ETOPO05_X'
        stry = 'ETOPO05_Y'
        strb = 'ROSE'
        #FIXME
        #This only works for ocidental hemisphere
        blon = (360-etopodata.variables[strx][:])*-1
        blat = etopodata.variables[stry][:]
    else:
        raise ValueError('Selected etopo not configured yet')
        

    
    condy = (blat>ltd)&(blat<ltu)
    condx = (blon>lnd)&(blon<lnu)
    
    topoin = etopodata.variables[strb][condy,:][:,condx]
    
    bLON,bLAT = np.meshgrid(blon[condx],blat[condy])
    
    if just_bathy:
        return bLON,bLAT,topoin
        
    else:
        plt.ioff()
        plt.figure()
        C = plt.contour(bLON,bLAT,-topoin,levels=[-PROF],colors='k',linewidths=1)
        plt.close()
        C = np.vstack(C.allsegs[0])
        pseudo_lat = C[:,1]
        pseudo_lon = C[:,0]
        pseudo_dat = np.zeros(pseudo_lat.shape)
        return pseudo_lon[::step],pseudo_lat[::step],pseudo_dat[::step]            

    

