#This program reads .bl files and .cnv files and extract it
#to xls table

#importing packages
from glob import glob
import numpy as np
import pandas as pd
import io

#running extra functions
%run /data0/iury_data/Copy/ocean/CTD/ctdproc.py

#defining functions
def scan_btl(fname,skiprows=2):
    names = ['bottle_order','bottle_no','time',
                'scan_i','scan_f']
    dat = pd.read_csv(fname,skiprows=skiprows,names=names,delimiter=',')
    
    return dat['bottle_order'].values,dat['scan_f'].values
    
def ctd_lonlat_bottlesprof(fname,sta,listactd):
    btl,scan_f = scan_btl(fname,skiprows=2)
    
    fctd = listactd[sta==int(fname.split('_')[2])]
    
    lon,lat,names,skiprows = get_lonlatnames(fctd[0],
                            lathint='Latitude =',
                            lonhint='Longitude =')
    CTD = np.array(pd.read_csv(fctd[0],skiprows=skiprows,
 			names=names,delim_whitespace=True,
 			usecols=['scan','prDM']))
    
    profbl = []
    for scan in scan_f:
        profbl.append(CTD[CTD[:,0]==scan,1])
    profbl = np.hstack(profbl)

    line = np.zeros(27)*np.nan
    line[0] = int(fname.split('_')[2])
    line[1] = lon
    line[2] = lat
    line[(btl+2).tolist()] = profbl
    return line
    
#Setting paths
                        
#path = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL02/BL'
#pathctd = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL02/CTD_hann'
#pathsave = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL02/'

pathsave = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/OMP/'
path = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/BL'
pathctd = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/CTD_hann'
pathsave = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/'

#reading .bl files
lista = glob(path+'/*.bl')
#sorting it
lista.sort()

#reading .cnv files
listactd = glob(pathctd+'/*.cnv')
#sorting it
listactd.sort()
#converting list to array
listactd = np.array(listactd)

#geting the station number
sta = []
for f in listactd:
    sta.append(int(f.split('_')[3]))
sta = np.array(sta)

#creating table numeric data
tab = []
for fname in lista:
    print fname.split('/')[-1].split('.')[0]
    try:
        tab.append(ctd_lonlat_bottlesprof(fname,sta,listactd))
    except:
        print 'No bottles on '+fname.split('/')[-1].split('.')[0]
TAB = np.vstack(tab)

#Creting table string columns data
names = ['Station','Longitude','Latitude']
for i in np.arange(24)+1:
    names.append('Bottle%.2i'%(i))
    
#Creating DataFrame
DF = pd.DataFrame(TAB,columns=names)

#saving the DataFrame
DF.to_excel(pathsave+'bottle.xls','io.excel.xls.writer')





