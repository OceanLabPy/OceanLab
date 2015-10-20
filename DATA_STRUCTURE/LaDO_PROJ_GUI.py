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
import time as tm
import getpass

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
    root.wm_geometry(newGeometry='1100x150+500+500')
    
    Label(root, text="New Project?").grid(sticky='w',row=0)
    Label(root, text="Add /NAMEOFPROJECT in the end").grid(sticky='w',row=1)
    Label(root, text="Path to Project Directory").grid(sticky='w',row=2)
    
    
    path = Entry(root,width=30)
    path.insert(END,os.getcwd())
    path.grid(row=3,column=0,sticky='w')
    
    options = ["CRUISE","GLIDER","MODELLING","MOORING"]
        
    v = IntVar()
    v.set(1)
    cnt = 4
    for txt in options:
        button_str = "Radiobutton(root,text=txt,command=v.set(cnt),variable=v,value=cnt)"
        exec('button = '+button_str)
        button.grid(row=2,column=cnt)
        cnt+=1
 
    Label(root, text="Period (MONTH/YEAR)").grid(sticky='w',row=16)
    date_i = Entry(root,width=20)
    date_i.grid(row=17,column=0,sticky='w')
    Label(root, text="to").grid(row=17,column=1)
    date_f = Entry(root,width=20)
    date_f.grid(row=17,column=2,sticky='w')
              
    Button(root, text='XINHO', command=root.quit).grid(row=18, sticky=W, pady=4)
    mainloop( )
    date_i,date_f,path_out,option = date_i.get(),date_f.get(),path.get(),v.get()
    root.destroy()
    
    return date_i,date_f,path_out,option


# Defines the path of the folder project and the properties that will be added
date_i,date_f,path,v = GUI()
tm=tm.localtime()
d,m,y,h,mn,s=tm.tm_mday,tm.tm_mon,tm.tm_year,tm.tm_hour,tm.tm_min,tm.tm_sec
user,machine=getpass.getuser(),os.uname()[1]

# Tests if the user input the / at the end of the directory
if path[-1]=='/':
    path = path[:-1]
    
# Changes the path in order to make sure that all the folders are in uppercase    
path=os.path.dirname(path)+'/'+path.split('/')[-1].upper()

#Creates the .pdf folder with the tree system
seeme1 = os.path.join(script_dir,'SEEME.cmap')
seeme2 = os.path.join(script_dir,'SEEME.pdf')


try:
    # checks if the folder already exists. if the project folder already exists the fuction only updates.
    os.chdir(path)
    print 'Project Directory Updated!'
except:
    # if not creates a new project folder
    print 'Project Directory Created!'
    os.system('mkdir '+path)
    os.system('mkdir '+os.path.join(path,'PLAN'))
    os.system('mkdir '+os.path.join(path,'PROJECT_REPORT'))
    #Creatig the README.txt file
    f=open(os.path.join(path,'README.txt'),'w')
    safe=open(os.path.join(path,'.README_s.txt'),'w')
    text='''

Created in %s/%s/%s at %s:%s:%s by: %s@%s

This is a folder-tree pattern created and used by Ocean Dynamic Laboratory from the Oceanographic Institute - University of São Paulo.
The scheme of the tree system used to create the folders is in the SEEME.pdf file.

In this folder system there are all the data and processing routines used in the %s project.

'''%(d,m,y,h,mn,s,user,machine,path.split('/')[-1].upper(),)
    f.write(text)
    safe.write(text)
    f.close()
    safe.close()
    os.system('cp '+seeme1+' '+path)
    os.system('cp '+seeme2+' '+path)  

#date_i='10/2015' ##### TODO
#date_f='11/2015' ##### TODO

if v==1:
    name,datanames = CRUISE_GUI()
    os.system('mkdir '+os.path.join(path,name))
    f=open(os.path.join(path,'README.txt'),'a')
    safe=open(os.path.join(path,'.README_s.txt'),'a')
    text='''

Edited in %s/%s/%s at %s:%s:%s by: %s@%s

In this project there was a cruise from %s to %s named %s. 
The parameter collected in the cruise were:'''%(d,m,y,h,mn,s,user,machine,date_i,date_f,name)
    f.write(text)
    safe.write(text)
    for prop in datanames:
        txt='''
        -%s'''%(prop)
        f.write(txt)
        safe.write(txt)
    f.close()
    safe.close()
    fig_dirs(os.path.join(path,name))
    print name
    for datname in datanames:
        print datname
        dataproc_dirs(datname,os.path.join(path,name))
elif v==2:
    f=open(os.path.join(path,'README.txt'),'a')
    safe=open(os.path.join(path,'.README_s.txt'),'a')
    text='''

Edited in %s/%s/%s at %s:%s:%s by: %s@%s

In this project it was done a glider sampling from %s to %s. '''%(d,m,y,h,mn,s,user,machine,date_i,date_f)
    f.write(text)
    safe.write(text)
    f.close()
    safe.close()
    dataprocfigrep_dirs('GLIDER',path,reportout=False)
elif v==3:
    os.system('mkdir '+os.path.join(path,'MODELLING'))
    f=open(os.path.join(path,'README.txt'),'a')
    safe=open(os.path.join(path,'.README_s.txt'),'a')
    text='''

Edited in %s/%s/%s at %s:%s:%s by: %s@%s
    
In this project it was done a modelling experiment with output data from %s to %s. '''%(d,m,y,h,mn,s,user,machine,date_i,date_f)
    f.write(text)
    safe.write(text)
    f.close()
    safe.close()
    dataprocfigrep_dirs('OUTPUT',os.path.join(path,'MODELLING'),raw_proc=False)
elif v==4:
    f=open(os.path.join(path,'README.txt'),'a')
    safe=open(os.path.join(path,'.README_s.txt'),'a')
    text='''

Edited in %s/%s/%s at %s:%s:%s by: %s@%s
    
In this project it was done a mooring sampling from %s to %s. '''%(d,m,y,h,mn,s,user,machine,date_i,date_f)
    f.write(text)
    safe.write(text)
    f.close()
    safe.close()
    dataprocfigrep_dirs('MOORING',path,reportout=False)

