# -*- coding: UTF-8 -*-
#THIS FUNCTION CREATES THE DIRECTORY STRUCTURE OF DEFINED PROJECT
#IN THE CURRENT DIRECTORY. THE STRUCTURE AGREE WITH THE FLUXOGRAM
#CALLED SEEME.pdf/.cmap/.svg/.ps AND WAS MADE TO BE THE PATTERN
#OF OCEANOGRAPHIC PROJECT PROCESSED BY LaDO (Laboratório de Dinâmica
#Oceânica) OF UNIVERSITY OF SÃO PAULO.

#STRUCT MADE BY: Iury Sousa, Dante Campagnoli and Hélio Almeida
#SCRIPT MADE BY: Iury Sousa (simoesiury@gmail.com)

import os
from glob import glob




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


def dataprocfigrep_dirs(name,path,reportout=True):
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
    os.system('mkdir '+os.path.join(path,name,'DATA','RAW'))
    os.system('mkdir '+os.path.join(path,name,'DATA','PROCESSED'))
    
    os.system('mkdir '+os.path.join(path,name,'PROC'))
    
        

       
try:
    project_dir = input('Whats the project path?   ')
except:
    raise ValueError('ERROR: Did you forget to put \" or \' ??')
    
if project_dir[-1]!='/':
    project_dir += '/'


try:
    project_name = input('Whats the project name?   ')
except:
    raise ValueError('ERROR: Did you forget to put \" or \' ??')

seeme1 = os.path.join(script_dir,'SEEME.cmap')
seeme2 = os.path.join(script_dir,'SEEME.pdf')

os.system('mkdir '+os.path.join(project_dir,project_name))
os.system('mkdir '+os.path.join(project_dir,project_name,'PLAN'))
os.system('mkdir '+os.path.join(project_dir,project_name,'PROJECT_REPORT'))
        
os.system('cp '+seeme1+' '+os.path.join(project_dir,project_name))
os.system('cp '+seeme2+' '+os.path.join(project_dir,project_name))  

os.chdir(os.path.join(project_dir,project_name))

expcond = input('There are any expeditions?   True/False    ')
if expcond:
    expeditions = []
    
    cnt= 1
    while True:
        try:
            expeditions.append(input('Whats the name of expedition %.i?   '%(cnt)))
        except:
            raise ValueError('ERROR: Did you forget to put \" or \' ??')
            
        if input('There are more expeditions?   True/False    '):
            None
        else:
            break
        cnt += 1
 
     
expdataname = ['XBT','SAT','WATER','LADCP','ADCP','TERMOSAL']
if expcond:
    for exp in expeditions:
        path = os.path.join(project_dir,project_name,exp)
        os.system('mkdir '+path)

        fig_dirs(path)
        
        #os.chdir(path)
        
        for datname in expdataname:
            if input('There are any %s data?   True/False    '%(datname)):
               dataproc_dirs(datname,path) 
        
        #os.system('mkdir '+os.path.join(project_dir,project_name))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    