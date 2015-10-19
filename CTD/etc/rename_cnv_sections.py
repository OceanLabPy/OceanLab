import os


path = '/data0/iury_data/Copy/Dados/ONE1/CTD_hann/'
subfolders = [x[0] for x in os.walk(path)][1:]
subfolders.sort()

for folder in subfolders:
    os.chdir(folder)
    radname = folder.split('/')[-1]
    os.system('rename s/\'.cnv\'/\'_%s.cnv\'/g *'%(radname))