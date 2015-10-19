# -*- coding: UTF-8 -*-
#THIS FUNCTION CREATES THE DIRECTORY STRUCTURE OF DEFINED PROJECT
#IN THE CURRENT DIRECTORY. THE STRUCTURE AGREE WITH THE FLUXOGRAM
#CALLED SEEME.pdf/.cmap/.svg/.ps AND WAS MADE TO BE THE PATTERN
#OF OCEANOGRAPHIC PROJECT PROCESSED BY LaDO (Laboratório de Dinâmica
#Oceânica) OF UNIVERSITY OF SÃO PAULO.

#STRUCT MADE BY: Iury Sousa, Dante Campagnoli and Hélio Almeida
#SCRIPT MADE BY: Iury Sousa (simoesiury@gmail.com) and Hélio Almeida

import os
from glob import glob
import tkFileDialog as filedialog
from Tkinter import *

script_dir = os.path.dirname(os.path.realpath(__file__))


def dataproc_dirs(name,path):
    os.system('mkdir '+os.path.join(path,name))
    
    os.system('mkdir '+os.path.join(path,name,'DATA'))
    os.system('mkdir '+os.path.join(path,name,'DATA','RAW'))
    os.system('mkdir '+os.path.join(path,name,'DATA','PROCESSED'))
    
    os.system('mkdir '+os.path.join(path,name,'PROC'))

def fig_dirs(path):
    os.system('mkdir '+os.path.join(path,'FIG'))
    
    os.system('mkdir '+os.path.join(path,'FIG','SECTION'))
    os.system('mkdir '+os.path.join(path,'FIG','MAP'))
    os.system('mkdir '+os.path.join(path,'FIG','PROFILE'))
    os.system('mkdir '+os.path.join(path,'FIG','OTHER'))


def dataprocfigrep_dirs(name,path,reportout=True,raw_proc=True):
    os.system('mkdir '+os.path.join(path,name))
    
    if reportout:
        os.system('mkdir '+os.path.join(path,'REPORT'))    
    else:
        os.system('mkdir '+os.path.join(path,name,'REPORT'))

    os.system('mkdir '+os.path.join(path,name,'FIG'))
    
    os.system('mkdir '+os.path.join(path,name,'FIG','SECTION'))
    os.system('mkdir '+os.path.join(path,name,'FIG','MAP'))
    os.system('mkdir '+os.path.join(path,name,'FIG','PROFILE'))
    os.system('mkdir '+os.path.join(path,name,'FIG','OTHER'))
    
            
            
    os.system('mkdir '+os.path.join(path,name,'DATA'))
    if raw_proc:
        os.system('mkdir '+os.path.join(path,name,'DATA','RAW'))
        os.system('mkdir '+os.path.join(path,name,'DATA','PROCESSED'))
        
    os.system('mkdir '+os.path.join(path,name,'PROC'))
    
def CRUISE_GUI():    
    frame = Tk()
    
    frame.wm_title('Cruise')
    frame.wm_geometry(newGeometry='110x250+500+500')
    

    Label(frame, text="Cruise name:").grid(sticky='w',row=0)
    name = Entry(frame,width=10)
    name.grid(row=1,column=0,sticky='w')
        
    
    options = ["ADCP","CTD","LADCP","SAT","XBT","WATER"]
    
    cnt = 2
    for txt in options:
        exec('%s = IntVar()'%(txt))
        
        button_str = "Checkbutton(frame,onvalue=cnt,text=txt,variable=%s)"%(txt)
        exec('button = '+button_str+'.grid(row=cnt,column=0,sticky=\'w\',pady = 4)')
        cnt+=1
    
    Button(frame, text='XINHO', command=frame.quit).grid(row=18, sticky=W, pady=4)
    mainloop( )
    
    opts = []
    for op in options:
        exec('cond = %s.get()!=0'%(op))
        if cond:
            opts.append(op)
            
    return name.get().upper(),opts


def GUI():    
    root = Tk()
    root.wm_title('PROJECT')
    root.wm_geometry(newGeometry='600x110+500+500')
    
    Label(root, text="Path to Project Directory").grid(sticky='w',row=0)
    path = Entry(root,width=20)
    path.grid(row=1,column=0,sticky='w')
    
    options = ["CRUISE","GLIDER","MODELLING","MOORING"]
        
    v = IntVar()
    v.set(1)
    cnt = 1
    for txt in options:
        button_str = "Radiobutton(root,text=txt,command=v.set(cnt),variable=v,value=cnt)"
        exec('button = '+button_str)
        button.grid(row=2,column=cnt)
        cnt+=1
    
    Button(root, text='XINHO', command=root.quit).grid(row=18, sticky=W, pady=4)
    mainloop( )
    path_out,option = path.get(),v.get()
    root.destroy()
    
    return path_out,option

path,v = GUI()

seeme1 = os.path.join(script_dir,'SEEME.cmap')
seeme2 = os.path.join(script_dir,'SEEME.pdf')

try:
    os.chdir(path)
    print 'Project Directory Updated!'
except:
    print 'Project Directory Created!'
    os.system('mkdir '+path)
    os.system('mkdir '+os.path.join(path,'PLAN'))
    os.system('mkdir '+os.path.join(path,'PROJECT_REPORT'))
    
    os.system('cp '+seeme1+' '+path)
    os.system('cp '+seeme2+' '+path)  


if v==1:
    name,datanames = CRUISE_GUI()
    
    os.system('mkdir '+os.path.join(path,name))
    fig_dirs(os.path.join(path,name))
    print name
    for datname in datanames:
        print datname
        dataproc_dirs(datname,os.path.join(path,name))
elif v==2:
     dataprocfigrep_dirs('GLIDER',path,reportout=False)
elif v==3:
     os.system('mkdir '+os.path.join(path,'MODELLING'))
     dataprocfigrep_dirs('OUTPUT',os.path.join(path,'MODELLING'),raw_proc=False)
elif v==4:
     dataprocfigrep_dirs('MOORING',path,reportout=False)

