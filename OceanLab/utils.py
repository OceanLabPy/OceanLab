import numpy as np

##### Function
#=============================================================================
# NEAREST DISTANCE
#=============================================================================
def argdistnear(x,y,xi,yi):
    '''
    This function finds the index to nearest points in (xi,yi) from (x,y).
    
    usage:
    x,y = [5,1,10],[2,6,3]
    xi,yi = np.linspace(0,19,20),np.linspace(-5,30,20)
    ind = argdistnear(x,y,xi,yi)

    INPUT:
    --> (x,y): points [list]
    --> (xi,yi): series to search nearest point [list]
    '''
    idxs = [np.argmin(np.sqrt((xi-xx)**2 + (yi-yy)**2)) for xx,yy in zip(x,y)]
    idxs = np.array(idxs)
    return idxs
#=============================================================================
