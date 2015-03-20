from ctd import DataFrame, movingaverage
import os
from ocean import *


#Using the package ocean and ctd you can export the processed CTD data
#to pickle format that can be joint as a section by the function return_section()
#from ocean.seaext
#Part of this script is based on orientations given in:
#    https://ocefpaf.github.io/python4oceanographers/blog/2013/07/29/python-ctd/

#define the args
kw = dict(below_water=True,lon=13,lat=14)
#file path
fname = '/home/iury/ocean/Examples/CTD_data_example.cnv'
#reading cnv
cast = DataFrame.from_cnv(fname, **kw)
#split downcast and upcast
downcast, upcast = cast.split()
#getting the longitude and latitude
lon,lat = get_lonlat(fname)
#using just downcast
data = downcast

#very basic processing with window based on data size
if data.shape[0]<101: # se o tamanho do perfil for com menos de 101 medidas

	if (data.shape[0]/2)%2 == 0: # caso a metade dos dados seja par
		blk = (data.shape[0]/2)+1 # bloco = a metade +1
	else:
		blk = data.shape[0]/2 # se for impar o bloco e a metade
		# remove spikes dos perfis de temperatura e condutividade
	temp = downcast['t068C'].despike(n1=2, n2=20, block=blk)

	temp = temp.bindata(delta=1)
	temp = temp.interpolate()

	condt = downcast['c0S/m'].despike(n1=2, n2=20, block=blk)

	condt = condt.bindata(delta=1)
	condt = condt.interpolate()
		
else:
	blk=101
	# para perfis com mais de 101 medidas, utiliza-se blocos de 101
	temp = downcast['t068C'].despike(n1=2, n2=20, block=blk)

	temp = temp.bindata(delta=1)
	temp = temp.interpolate()

	condt = downcast['c0S/m'].despike(n1=2, n2=20, block=blk)

	condt = condt.bindata(delta=1)
	condt = condt.interpolate()

#convert to Dataframe
df = pd.DataFrame(index=temp.index)
#calculating the new variables
df['t090'] = gsw.t90_from_t68(temp)
df['sp'] = gsw.SP_from_C(condt*10,temp,df.index.values)
df['sa'] = gsw.SA_from_SP(df['sp'],df.index,lon,lat)
df['pt'] = gsw.pt_from_t(df['sa'],temp,df.index.values)
df['ct'] = gsw.CT_from_pt(df['sa'],df['pt'])
df['dens0'] = sw.dens0(df['sp'],temp)
df['psigma0'] = sw.pden(df['sp'],temp,df.index.values,pr=0)-1000
df['psigma1'] = sw.pden(df['sp'],temp,df.index.values,pr=1000)-1000
df['gpan'] = sw.gpan(df['sp'],temp,df.index.values)
df['lat'] = lat
df['lon'] = lon

#Now you export as a pickle the data
pd.to_pickle(df,fname.split('.')[0])

#In the END you can use the function return_section to read all
#the stations and join as a section