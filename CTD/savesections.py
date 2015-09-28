import pandas as pd
import scipy.io as sio
from glob import glob

def savesection(path,ext):
    section = return_section(path,'*'+ext+'.cnv')
    #salvando pickle
    section.to_pickle(path+'sections/'+ext)
    ##salvando excel
    #section.to_excel(path+'all/alldata.xls','io.excel.xls.writer')
    #salvando em matlab
    Matlab = dict()
    for prop in section.minor_axis:
        Matlab.update({prop:section.minor_xs(prop).values})
    Matlab.update({'PRESS':section.major_axis.values})
    sio.savemat(path+'sections/'+ext,Matlab)


#LENDO TUDO
path = '/home/iury/Copy/Dados/OL2/CTD_hann/'
section = return_section(path,'*.cnv')
#salvando pickle
section.to_pickle(path+'all/alldata')
##salvando excel
#section.to_excel(path+'all/alldata.xls','io.excel.xls.writer')
#salvando em matlab
Matlab = dict()
for prop in section.minor_axis:
    Matlab.update({prop:section.minor_xs(prop).values})
Matlab.update({'PRESS':section.major_axis.values})
sio.savemat(path+'all/alldata',Matlab)


#LENDO TUDO
lista = glob(path+'*.cnv')

rads = []
for i in lista:
    rads.append(i.split('/')[-1].split('.')[0].split('_')[-1])
rads = np.unique(np.array(rads))

for rad in rads:
    savesection(path,rad)







