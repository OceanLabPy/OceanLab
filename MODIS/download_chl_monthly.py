from matplotlib.dates import date2num
from datetime import datetime
from calendar import monthrange
import numpy as np
import os


def year_month_lims(Y):
    inis,fins = [],[]
    cnt = 1
    for mo in np.arange(12)+1:
        mrange = monthrange(Y,mo)[1]
        if mo==1:
            inis.append(1)
        else:
            inis.append(np.sum(np.subtract(fins,inis))+cnt)
        
        fins.append(inis[cnt-1]+mrange-1)
        
        cnt+=1
    
    lims = zip(inis,fins)
    return lims

def download_chl_monthly(Y,path,datatype='A',res='9km'):
    os.chdir(path)
    
    errlist = []
    
    lims = year_month_lims(Y)
    
    baseurl = 'http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/'
    suffix = '.L3m_MO_CHL_chlor_a_'+res+'.nc'
    
    month = 1
    for ini,fin in lims:
        basefile = datatype+'%.4i%.3i%.4i%.3i'%(Y,ini,Y,fin)
        proc = os.system('wget '+baseurl+basefile+suffix)
        if proc == 0:
            print 'Month %.2i of %.4i Downloaded!'%(month,Y)
        else:
            errlist.append(baseurl+basefile+suffix)
        month += 1
        
    return errlist
        
path = '/data0/iury_data/Copy/Dados/MODIS/CHL/A/'


errall = []

for yr in np.arange(2002,2016):
    errall.append(download_chl_monthly(yr,path,res='4km'))   

for err_urls in errall:
    if err_urls!=[]:
        for err_url in err_urls:
            proc = os.system('wget '+err_url)
        
        
        
