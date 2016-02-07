#This program reads the data from chemical analysis
#by Chemical Oceanography Laboratory from USP
#and the data extracted by bottles.py and
#creates a table with chemical analysis, longitude,
#latitude and depth of amostration

#importing packages
import pandas as pd
import numpy as np
from collections import Counter,OrderedDict


def return_tempsal(ctdfile,dpth):
    CTD = pd.read_pickle(ctdfile[:-4])
    ox = (CTD['sbeox0Mm/Kg']/((CTD['psigma0']+1000)/1000))
    ox = ox[dpth].values
    return CTD['pt'][dpth].values,CTD['sp'][dpth].values,ox

#seting paths
expno=1
pathchem = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/OMP/data/'
pathbl = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/OMP/data/'
pathsave = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/OMP/data/nutrients/'
pathctd = '/data0/iury_data/Copy/Dados/MARSEAL/MARSEAL01/CTD_hann'


#reading .cnv files
listactd = glob(pathctd+'/*.cnv')
#sorting it
listactd.sort()
#converting list to array
listactd = np.array(listactd)

#geting the station number
stactd = []
for f in listactd:
    stactd.append(int(f.split('_')[3]))
stactd = np.array(stactd)

#define the masses number
Dm = OrderedDict({'AT':0,'ACAS':1,'AIA':2,'ACS':3,'APANS':4,
'APANI':5,'ACI':6,'AFA':7,'AS':8,'EXTRA':9})

index = ['Longitude','Latitude','Depth',
         'Nitrito','NitratoNitrito','Fosfato',
         'Silicato','TemperaturaPotencial','Salinidade','Oxigenio']
columns = ['AT','ACAS','AIA','ACS','APANS',
'APANI','ACI','AFA','AS','EXTRA']

#Readind data
dat = pd.read_csv(pathchem+'MARSEAL_I_PY.csv')
datbl = pd.read_csv(pathbl+'bottle.csv')

#Getting the station numbers and water masses analysed
sta = []
wm = []
for name in dat['Amostra']:
    sta.append(int(name.split(' - ')[0][1:]))
    wm.append(name.split(' - ')[1])
sta,wm = np.array(sta),np.array(wm)


#Counting the number of masses analysed by each station
C = np.array(Counter(sta).items())

#Geting the .bl data about depth of amostration
profs = datbl.iloc[:,4:].values

#UFF nan
profs[profs<=10] = np.nan

#Geting the .bl data about depth of stations
stabl = datbl['Station'].values




#creating table numeric data
DATA = OrderedDict()
for st,masses in zip(C[:,0],C[:,1]):
    
    wmpos = [Dm[wms] for wms in wm[np.array(sta)==st].tolist()]
    
    tab = np.zeros((10,10))*np.nan
    
    p = profs[stabl==st,:][0]
    difp = np.abs(np.diff(p))
    difp[np.isnan(difp)] = 0
    indx = np.argsort(difp)
    pmean = []
    lims = np.hstack([-1,np.sort(indx[-(masses-1):]),p.size])
    for i,f in zip(lims[:-1],lims[1:]):
        pmean.append(np.nanmean(p[i+1:f+1]))
    pmean = pmean[::-1]
    
    tab[2,wmpos] = np.round(pmean)
    tab[0,wmpos] = datbl['Longitude'][stabl==st]
    tab[1,wmpos] = datbl['Latitude'][stabl==st]
    tab[3:7,wmpos] = dat.iloc[sta==st,1:].values.T
    
    temp,sal,oxi = return_tempsal(listactd[stactd==st][0],
                                np.round(pmean))
    
    tab[7,wmpos] = temp
    tab[8,wmpos] = sal
    tab[9,wmpos] = oxi
    
    DF = pd.DataFrame(tab,columns=columns,index=index)
    
    DF.to_csv(pathsave+'SEAL%.2i_%.2i'%(expno,st))
    
    DATA.update({'SEAL%.2i_%.2i'%(expno,st):DF})

DATA = pd.Panel.from_dict(DATA)
    


    
    
    
    
    
    
    
    
