# -*- coding: utf-8 -*-
import numpy as np
import seawater as sw


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
