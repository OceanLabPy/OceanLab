# -*- coding: utf-8 -*-
#These functions there's no general use. I maintain and publish this for some reason that I donno
#Made by Iury Sousa - Sao Paulo/Brazil


def angcalc(x,y):
	rads = []
	dx,dy = np.diff(x),np.diff(y)
	for i,j in zip(dy,dx):
		rads.append(atan2(i,j))
	return rads


def quadrant(angl,degrees=False):
	
	oang = angl
	
	if degrees:
		angl = np.radians(angl)			
	while angl<0:
		angl+= 2*np.pi	
	while angl>2*np.pi:
		angl -= 2*np.pi	
	if (angl>0)&(angl<np.pi/2.):
		qrt = 1	
	elif (angl>np.pi/2)&(angl<np.pi):
		qrt = 2	
	elif (angl>np.pi)&(angl<3*np.pi/2.):
		qrt = 3	
	elif (angl>3*np.pi/2.)&(angl<2*np.pi):
		qrt = 4	
	return qrt



def westraxll(exped):
#	filepath = '/media/lof/SAMSUNG/tcc/dados/WESTRAX/WESTRAX%i/PEG' %(exped)

	filepath = '/home/iury/Copy/TCC/dados/WESTRAX/WESTRAX%i/PEG' %(exped)

	profiles = np.sort(glob(os.path.join(filepath,'*.pro')))

	lats = []
	lons = []
	prof = []

	U = OrderedDict()
	V = OrderedDict()
	for pro in profiles:
		head,exp,est,lat,lon,upro,vpro,exdata = [],[],[],[],[],[],[],[]
		data = np.loadtxt(pro)
		head = data[0,:]

		if np.all(data[1::,7]==-9999):
			print('Estação %s só contém NaNs'%(pro.split('/')[-1].split('.')[0]))	
		else:	
			lat = data[0,4]
			lon = -data[0,5]
	
			mask = data[1::,7]!=-9999
			u = data[1::,7][mask]/100
			v = data[1::,8][mask]/100
			
			uinterp = scint.interp1d(data[1::,1][mask],u)
			vinterp = scint.interp1d(data[1::,1][mask],v)
			
			u1 = uinterp(np.arange(0,data[1::,1][mask].max()))		
			v1 = vinterp(np.arange(0,data[1::,1][mask].max()))
			
			nna = np.tile(np.nan,5050-u1.size)
			u1 = np.hstack((u1,nna))
			v1 = np.hstack((v1,nna))	


			U.update({pro.split('/')[-1].split('.')[0]: u1})
			V.update({pro.split('/')[-1].split('.')[0]: v1})
			
			lats.append(lat)
			lons.append(lon)
		
	U = pd.DataFrame.from_dict(U)
	V = pd.DataFrame.from_dict(V)
	
	lons = np.array(lons)
	lats = np.array(lats)
	
	return lons,lats,U.values,V.values
	

#################################################################

def adcpnewll(pathfile):
	data = sio.loadmat(pathfile)
	
	mask = ~np.isnan(data['lat'][0])
	lats = data['lat'][0][mask]
	lons = data['lon'][0][mask]-360
	u = data['u'][:,mask]
	v = data['v'][:,mask]
	prf = data['prf']
	
	dist = np.hstack((0,np.cumsum(sw.dist(lats,lons)[0])))

	u = pd.DataFrame(u,index=prf,columns=dist)
	v = pd.DataFrame(v,index=prf,columns=dist)
	
	u2 = interpgrid(u,interp=False,dya=600)
	v2 = interpgrid(v,interp=False,dya=600)
		
	return lons,lats,u2.values,v2.values
	