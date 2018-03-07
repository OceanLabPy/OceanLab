# -*- coding: utf-8 -*-
from scipy.optimize import curve_fit
import itertools
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from OceanLab.utils import *

def vectoa(Xg,Yg,X,Y,U,V,corrlenx,corrleny,err,b=0):
    '''
    Vectoa is a vectorial objective analysis function.

    It interpolates a velocity field (U and V, east and north velocity components)
    into a streamfunction field (dimension MxN) solving the Laplace equation:

        $nabla^{2}Psi=0$.

   ======================================================================

    Input:
        Xg & Yg - Grid of of interpolation points (i.e. LON & LAT grid)

        X & Y - Arrays of observation points

        U & V - Arrays of observed east and north components of velocity

        corrlen & err - Correlation length scales (in x and y) and error for a
                    gaussian streamfunction covariance function (floats)

        b - Constant value that forces a correction in the data mean value.
            Unless it is defined, b=0 by default.


   ======================================================================
    Output:

        PSI - Streamfuction field matrix with MxN dimension.
              The dimension of the output is always the same of XC & YC

   ======================================================================

    PYTHON VERSION by:
       Iury Sousa and Hélio Almeida - 30 May 2016
       Laboratório de Dinâmica Oceânica - IOUSP

   ======================================================================'''
   # making sure that the input variables aren't changed
    xc,yc,x,y,u,v=Xg.copy(),Yg.copy(),X.copy(),Y.copy(),U.copy(),V.copy()

    corrlen = corrleny
    xc = xc*( corrleny*1./corrlenx)
    x = x*(corrleny*1./corrlenx)

    n = len(x)
    # Joins all the values of velocity (u then v) in one column-wise array.
    # Being u and v dimension len(u)=y, then uv.shape = (2*y,1)
    uv=np.array([np.hstack((u,v))]).T    #data entered row-wise


    # CALCULATING angles and distance

    # pp is a join of two matrix calculating the distance between every point of observation and all others like the
    # example below
    # len(y) = M -> pp[0].shape = MxM being:

    #    pp[0][i] = y-y[i]
    pp = -np.tile(y,(n,1)).T+np.tile(y,(n,1)),-np.tile(x,(n,1)).T+np.tile(x,(n,1))

    #
    t = []
    for ii,jj in zip(pp[0].ravel(),pp[1].ravel()):
        t.append(np.math.atan2(ii,jj))
    t = np.array(t)
    t.shape = pp[0].shape
    # t end up to be angles and d2 the distances between every observation point and all the others
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

    # angles and distances
    nv1,nv2 =xc.shape
    nv=nv1*nv2


    xc = xc.T.ravel()
    yc = yc.T.ravel()


    #% the same idea of pp but this time for interpolation points
    ppc = -np.tile(yc,(n,1)).T+np.tile(y,(nv,1)),-np.tile(xc,(n,1)).T+np.tile(x,(nv,1))
    tc = []
    for ii,jj in zip(ppc[0].ravel(),ppc[1].ravel()):
        tc.append(np.math.atan2(ii,jj))
    tc = np.array(tc)
    tc.shape = ppc[0].shape
    d2=((np.tile(xc,(n,1)).T-np.tile(x,(nv,1)))**2+(np.tile(yc,(n,1)).T-np.tile(y,(nv,1)))**2)
    R=np.exp(-lambd*d2)+bmo;


    P=np.zeros((nv,2*n))
    # streamfunction-velocity covariance
    P[:,0:n]=np.sin(tc)*np.sqrt(d2)*R;
    P[:,n:2*n]=-np.cos(tc)*np.sqrt(d2)*R;

    PSI=np.dot(P,np.linalg.solve(A,uv))   # solvi the linear system
    PSI=PSI.reshape(nv2,nv1).T

    return PSI

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
  xc = xc*( corrleny*1./corrlenx)
  x = x*(corrleny*1./corrlenx)

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

def fitting(cdat,dist,binn,cut=550):

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
    plt.ylabel(r'Correlation')

    ER = 1-func(0,*popt1)
    LC = popt1[-1]
    plt.text(LIM*0.5,0.4,
            r'$\epsilon^{2}$: %.3f - LC: %.3f km'%(ER,LC),
            fontsize=15)

    return ER,LC
