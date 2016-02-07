def my_eof_interp(M,nmodes,errmin=1e-15,repmax=None):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vi = my_eof_interp(M,nmodes)
    #%   This function try to fill gappy data using its EOFs
    #%
    #%   M is the data Matrix (generaly the diferent stations
    #%     	are on the columns);
    #%
    #%	nmodes is the number of statiscally significant EOF modes (to be used for the series reconstruction)  
    #%
    #%% The method -- based on Beckers & Rixen (2003, JAOTech)
    #%
    #% The algorithm proposed by the authors is as follows.  (see Eqs. (13) and (14) from Beckers & Rixen)
    #%
    #% (I)   Put zero in the GAPPY positions ("unbiased INITIAL GUESS, by assuming that the removed mean was unbiased). 
    #%
    #% (II)  Then you can have a first estimative of the EOFs.
    #%
    #% (III) You can reconstruct the series using the statiscally significant modes (e.g.: determined using Monte Carlo
    #%           simulations as in Preisendorfer [1988] -- OBS: the authors propose another method using cross validation 
    #%           during the iterative process)
    #% 
    #% (IV)  Substitute the reconstructed series ONLY in the GAPPY positions in your original matrix.
    #%
    #% (V)   Run steps (II), (III) and (IV) iteractively till its convergence (i.e. compare the reconstructed series at time
    #%           t+1 with that at time t).
    #%
    #% OBS: It doesn't work for filling gaps simultaneously in all depths.
    
    Mmean = np.nanmean(M,axis=1)
    
    M = M.T - np.nanmean(M,axis=1)
    M = M.T
    
    f = np.isnan(M)
    M[f] = 0 #%% initial guess in gappy positions
    
    rep = 0
    while True:
   	print rep+1 
    
    
   	Ud,Dd,Vd = np.linalg.svd(M,'econ')
   	Dd = np.diagflat(Dd)
   	
   	#% compute the explained variance (total variance is the sum of the Dd elements squared -total energy-)
   	evd = np.diag(Dd)**2 #%% remember that the singular values (Dd) are the sqrt of 
   	evd = Dd/np.sum(Dd)  #%%    the eigenvalue of the cov. matrix.
    
    
   	#%% truncated EOF reconstruction
   	Ud = Ud[:,0:nmodes]
   	Vd = Vd[:,0:nmodes]
   	Dd = Dd[0:nmodes,0:nmodes]
   	#Dd = Dd[0:nmodes]
    
   	#% reconstructing the series -- truncated reconstruction (SVD: Mrec = Ud * Dd *Vd')
   	Xa = np.dot(np.dot(Ud,Dd),Vd.T)
    
   	Mold = M.copy() #% saving for compute the error
        
        M[f] = Xa[f] 
    
   	#% testing the error
   	if rep>=1:
  		err = np.nanmean(np.abs(M - Mold))/(np.nanmax(M)-np.nanmin(M))
  		print err
  		if err < errmin:
 			break
        if repmax!=None:
            if rep>=repmax:
                break
        rep += 1

    #%% Saving the last the series with the EOF reconstructed series in the gappy positions;
    vi = (M.T+Mmean).T
    
    return vi