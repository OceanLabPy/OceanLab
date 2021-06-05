# -*- coding: utf-8 -*-
import numpy as np
import seawater as sw

# ===========================================================
# DYNAMICAL MODES AMPLITUDE
# ===========================================================
def vmode_amp(A,vi):
    '''
    This function makes the projection of every vertical mode to
    timeseries or section matrix to obtain its amplitude.
    It will be used to data reconstruction.

    Operation: ((A'*A)^-1)*(A'*vi)
    =============================================================
    
    INPUT:
       A  = 
       vi = 

    OUTPUT:
       AMP = 

    ==============================================================
    '''

    return np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,vi))

# ===========================================================
# RELATIVE VORTICITY
# ===========================================================
def zeta(x,y,U,V):
    '''
    ZETA k-component of rotational by velocity field
       Returns the scalar field of vorticity from U and V
       velocity components. X and Y are longitude and latitude
       matrices, respectively, as U and V. All matrices have
       same dimension. Grid could be curvilinear, but the error
       will increase.

       --> !!! IMPORTANT !!!  <-----------------------------------------
       The grid indexes IM and JM must have origin on the left-inferior
       corner, increasing to right and to up, respectively.

                                  GRID
                :            :            :            :
                |            |            |            |
       (JM = 2) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 1) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 0) o ---------- o ---------- o ---------- o --- ...
            (IM = 0)     (IM = 1)     (IM = 2)     (IM = 3)

    ================================================================
    USAGE:  ZETA = zeta(x,y,u,v)

    INPUT:
       x     = decimal degrees (+ve E, -ve W) [-180..+180]
       y     = decimal degrees (+ve N, -ve S) [- 90.. +90]
       u     = velocity zonal component [m s^-1]
       v     = velocity meridional component [m s^-1]

    OUTPUT:
       ZETA  = Relative vorticity field [s^-1]

    AUTHOR:
       Iury T. Simoes-Sousa and Wandrey Watanabe - 29 Jun 2016
       Laboratório de Dinâmica Oceânica - IOUSP
    ================================================================
    '''

    #função para cálculo de ângulos
    angcalc = lambda dy,dx: np.math.atan2(dy,dx)

    #funções para os diferenciais com borda
    dX      = lambda d: np.hstack([(d[:,1]-d[:,0])[:,None],
                        d[:,2:]-d[:,:-2],(d[:,-1]-d[:,-2])[:,None]])
    dY      = lambda d: np.vstack([(d[1,:]-d[0,:])[None,:],
                        d[2:,:]-d[:-2,:],(d[-1,:]-d[-2,:])[None,:]])

    # x = coordenada natural i
    # y = coordenada natural j
    # X = coordenada zonal cartesiana
    # Y = coordenada meridional cartesiana

    dyX = dX(y)*60*1852*np.cos(y*np.pi/180)
    dyY = dY(y)*60*1852

    dxX = dX(x)*60*1852*np.cos(y*np.pi/180)
    dxY = dY(x)*60*1852

    dsx = np.sqrt(dxX**2 + dyX**2)
    dsy = np.sqrt(dxY**2 + dyY**2)

    #Derivadas naturais de U
    dUy = dY(U)
    dUx = dX(U)

    #Derivadas naturais de V
    dVy = dY(V)
    dVx = dX(V)

    #a função map aplicará a função lambda angcalc em todos os elementos
    #indicados depois da virgula e armazenará em uma lista
    angX = list(map(angcalc,dyX.ravel(),dxX.ravel()))
    angY = list(map(angcalc,dyY.ravel(),dxY.ravel()))

    #fazendo rechape dos ângulos para o formato das matrizes de distâncias
    angX = np.reshape(angX,dsx.shape)
    angY = np.reshape(angY,dsy.shape)

    #calculando as componentes X e Y das derivadas
    dv1,dv2 =  (dVy/dsy)*np.cos(angY),  (dVx/dsx)*np.cos(angX)
    du1,du2 = -(dUy/dsy)*np.sin(angY), -(dUx/dsx)*np.sin(angX)


    # zerando valores aos quais a derivada tendeu ao infinito:
    # isto só acontece se uma das dimensões da grade for paralalela a um
    # eixo cartesiano
    dv1[np.isinf(dv1)] = 0
    dv2[np.isinf(dv2)] = 0
    du1[np.isinf(du1)] = 0
    du2[np.isinf(du2)] = 0

    #somando as componentes
    ZETA = dv1+dv2+du1+du2

    return ZETA

# ===========================================================
# U AND V FROM STREAMFUNCTION
# ===========================================================
def psi2uv(x,y,psi):
    '''
    PSI2UV Velocity components from streamfunction

      Returns the velocity components U and V
       from streamfunction PSI. X and Y are longitude and latitude 
       matrices, respectively, as PSI.
       All matrices have same dimension.

       --> !!! IMPORTANT !!!  <-----------------------------------
       The grid indexes IM and JM must have origin on the 
       lower-left corner, increasing to right and to up.

                                  GRID
                :            :            :            :
                |            |            |            |
       (JM = 2) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 1) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 0) o ---------- o ---------- o ---------- o --- ...
            (IM = 0)     (IM = 1)     (IM = 2)     (IM = 3)

    ==============================================================

    USAGE:  u,v = psi2uv(x,y,psi)

    INPUT:
       x    = decimal degrees (+ve E, -ve W) [-180..+180]
       y    = decimal degrees (+ve N, -ve S) [- 90.. +90]
       psi  = streamfunction [m^2 s^-1]

    OUTPUT:
       u    = velocity zonal component [m s^-1]
       v    = velocity meridional component [m s^-1]

    AUTHOR:
       Iury T. Simoes-Sousa and Wandrey Watanabe
       May 2016 - Laboratório de Dinâmica Oceânica - IOUSP
       Universidade de São Paulo
    ==========================================================
    '''

    #função para cálculo de ângulos
    angcalc = lambda dy,dx: np.math.atan2(dy,dx)

    #funções para os diferenciais com borda
    dX      = lambda d: np.hstack([(d[:,1]-d[:,0])[:,None],
                        d[:,2:]-d[:,:-2],(d[:,-1]-d[:,-2])[:,None]])
    dY      = lambda d: np.vstack([(d[1,:]-d[0,:])[None,:],
                        d[2:,:]-d[:-2,:],(d[-1,:]-d[-2,:])[None,:]])

    # x = coordenada natural i
    # y = coordenada natural j
    # X = coordenada zonal cartesiana
    # Y = coordenada meridional cartesiana

    dyX = dX(y)*60*1852*np.cos(y*np.pi/180)
    dyY = dY(y)*60*1852

    dxX = dX(x)*60*1852*np.cos(y*np.pi/180)
    dxY = dY(x)*60*1852

    dsX = np.sqrt(dxX**2 + dyX**2)
    dsY = np.sqrt(dxY**2 + dyY**2)

    dpx = dX(psi)
    dpy = dY(psi)

    #a função map aplicará a função lambda angcalc em todos os elementos
    #indicados depois da virgula e armazenará em uma lista
    angX = list(map(angcalc,dyX.ravel(),dxX.ravel()))
    angY = list(map(angcalc,dyY.ravel(),dxY.ravel()))
    #fazendo rechape dos ângulos para o formato das matrizes de distâncias
    angX = np.reshape(angX,dsX.shape)
    angY = np.reshape(angY,dsY.shape)

    #calculando as componentes u e v das velocidades calculadas em JM e IM
    v1,u1 =  (dpy/dsY)*np.cos(angY), -(dpy/dsY)*np.sin(angY)
    u2,v2 = -(dpx/dsX)*np.sin(angX),  (dpx/dsX)*np.cos(angX)

    # zerando valores aos quais a derivada tendeu ao infinito:
    # isto só acontece se uma das dimensões da grade for paralalela a um
    # eixo cartesiano
    v1[np.isinf(v1)] = 0
    v2[np.isinf(v2)] = 0
    u1[np.isinf(u1)] = 0
    u2[np.isinf(u2)] = 0

    #somando as componentes
    U = u1+u2
    V = v1+v2

    return U,V

# ===========================================================
# EQUATORIAL MODES
# ===========================================================
def eqmodes(N2,z,nm,lat,pmodes=False):
    '''
    This function computes the equatorial velocity modes

    ============================================================

    INPUT:
       N2     = Brunt-Vaisala frequency data array
       z      = Depth data array (equaly spaced)
       lat    = Latitude scalar
       nm     = Number of modes to be computed
       pmodes = If the return of pressure modes is required.
                 (Default is False)

    OUTPUT:
       Si     = Equatorial modes matrix with MxN dimension being:
                 M = z array size
                 N = nm
       Rdi    = Deformation Radii array
       Fi     = Pressure modes matrix with MxN dimension being:
                 M = z array size
                 N = nm
                 (Returned only if input pmodes=True)

    AUTHOR:
       Iury T. Simoes-Sousa, Hélio Almeida and Wandrey Watanabe
       2016 Laboratório de Dinâmica Oceânica
       Universidade de São Paulo
    ============================================================
    '''

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
    Si = np.array(list(map(onorm,Si.T))).T
    # to keep consistency multiply by plus/minus
    # to make the leading term positive
    Si = np.array(list(map(plus_minus,Si.T))).T

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
        Fi = np.array(list(map(onorm,Fi.T))).T
        # to keep consistency multiply by plus/minus
        # to make the leading term positive
        Fi = np.array(list(map(plus_minus,Fi.T))).T

        Si = np.hstack([sb,Si])
        Fi = np.hstack([fb,Fi])

        return Si,radii,Fi
    else:
        Si = np.hstack([sb,Si])
        return Si,radii


# ===========================================================
# DYNAMICAL MODES
# ===========================================================
def vmodes(N2,z,nm,lat,ubdy='N',lbdy='N'):

    """
    
    ========================================================
    INPUT:

    OUTPUT:

    AUTHOR:
    ========================================================
    """

    f0 = sw.f(lat)
    dz = np.abs(z[1] - z[0])

    # f2n2 array
    f2n2 = (f0**2)/N2

    # assembling matrices

    # N2F2 matrix
    N2F2 = np.eye(z.size-1)
    for i in range(0, z.size-1):
        N2F2[i,i] = (f2n2[i] + f2n2[i+1])/2

    # 1st order difference matrix
    A = np.eye(z.size)[:,:-1]
    for i in range(1,z.size):
        A[i,i-1] = -1	

    A = A/dz	

    # linear operator C = A*N2F2*(A.transpose())
    A = np.matrix(A)
    N2F2 = np.matrix(N2F2)

    C = A*N2F2*A.transpose()

    if ubdy == 'D':
        C = C[1:,1:]
    if lbdy == 'D':
        C = C[:-1,:-1]

    # solve the eigenvalue problem
    lam,fi = np.linalg.eig(C)

    i = np.argsort(lam)
    ei = lam[i]
    fi = fi[:,i]

    ei = ei[:nm]
    fi = fi[:,:nm]

    if ubdy == 'D':
        fi = np.append(np.matrix(np.zeros(nm)),fi, axis=0)
    if lbdy == 'D':
        fi = np.append(fi,np.matrix(np.zeros(nm)), axis=0)

    # normalizing to get orthonormal modes 
    H = np.abs(z).max()

    for i in range(nm):
        s = np.sqrt(dz*((np.array(fi[1:,i])**2 + np.array(fi[:-1,i])**2)/2/H).sum())
        fi[:,i] = fi[:,i]/s

    # to keep consistency multiply by plus/minus to make the leading term positive
    for i in range(nm):
        if np.sign(fi[1,i]) < 0:
            fi[:,i] = -fi[:,i]

    # compute the deformation radii [km]; the barotropic radius is computed by sqrt(gH)/f0
    radii = np.zeros(nm)
    radii[0] = np.sqrt(9.81*H)/np.abs(f0)/1000
    radii[1:] = 1/np.sqrt(ei[1:])/1000

    return fi, radii
