import matplotlib.cm as cm

lista = glob('/home/iury/Copy/Dados/OL2/CTD_hann/sections/*')
lista.sort()
lista = lista[::2]

colors = cm.rainbow(np.linspace(0, 1, len(lista)))

plt.figure()
for fname,color in zip(lista,colors):
    sec = pd.read_pickle(fname)
    lon = np.nanmean(sec.minor_xs('lon'),axis=0)
    lat = np.nanmean(sec.minor_xs('lat'),axis=0)
    for ln,lt,txt in zip(lon,lat,sec.items.values):
        plt.text(ln,lt,txt)
    plt.scatter(lon,lat,c=color,s=50,label=fname.split('/')[-1])

plt.legend(loc='best')
