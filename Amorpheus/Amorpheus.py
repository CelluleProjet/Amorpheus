# -*- coding: utf-8 -*-
"""

Amorpheus: a python-based software for the treatment of x-ray scattering data of amorphous and liquid systems
Version 1.0.5

Contacts: Silvia Boccato (silvia.boccato@cnrs.fr), Yiuri Garino (yiuri.garino@cnrs.fr)

Copyright (c) 2020-2022 Silvia Boccato, Yiuri Garino, Guillaume Morard and Chrystele Sanloup

Download: https://github.com/CelluleProjet/Amorpheus
How to cite: https://doi.org/10.1080/08957959.2022.2032032

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 724690).

"""

import sys

if hasattr(sys, 'ps1'):
    print("Running interactively")
    print()
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -sf')
    import sys
    import matplotlib.pyplot as plt
    terminal = False
else:
    import matplotlib.pyplot as plt
    print("Running in terminal")
    print()
    plt.ion() #activate interactive mode in terminal scripts
    terminal = True
import equations as eq
import numpy as np

import warnings
warnings.filterwarnings("ignore")
import time, traceback
from lmfit import Parameters, fit_report, minimize
from scipy import stats

from pathlib import Path
import os
from scipy.signal import find_peaks

"""
Input_File = True => Requires the file Init.txt
Input_File = False => Uses the values defined below
"""

Input_File = True 
Init_filename = 'Init.txt'
debug = False

"""
File Format

1 => Th No Background
2 => Th + Background
3 => Q No Background
4 => Q + Background

QStep ==> for debug purpose only in case of FileFormat == 2
0 => Real data from file
1 => constant dQ approximation 

S(Q) generation at low Q ==> to be implemented in future versions
0 => S0[0:diff] = 0
1 => S0[0:diff] = S0[diff]
2 => S0[0:diff] = lorentzian 

incoherent ==>
0 => background file was taken on sample
1 => background file was taken on empty cell

"""

FileFormat = 4
QStep = 0
SofQ = 1
incoherent = 1 

"""
Data Filename WITH EXTENSION
If no background use 'NoBk.xy'
"""
filename = 'Cerium.dat'
filename_bkg =  'Background.xy'

"""
Data analysis parameters
"""
bkfactor=  0.927
rhomi=34.68 #atoms/nm^3
rmin= 0.22 #Angstrom 
qmin = 15   
qmax = 78
lambd= 0.03738 #nm, used only in fileformat 1 or 2
Normalization = 0


""" Experiment parameters """
elements = ['Ce','Al','Ni','Cu','Fe','C']
ratio = np.array([0.70,0.10,0.10,.1,0,0])

    
""" Target Save Filename NO EXTENSION """
F3_save_filename = 'Norm'
F4_save_filename = 'MinRmin'
F5_save_filename = 'MinQmax'

""" Figures output resolution DPI (dots per inch) """
DPI = 300


"""

Here starts the program

"""


print('/\\'*32)
print("""
Welcome to Amorpheus: a python-based software for the treatment 
        of x-ray scattering data of amorphous and liquid systems 

Copyright(c) 2020-2022 Silvia Boccato, Yiuri Garino, 
                          Guillaume Morard and Chrystele Sanloup

Distributed under the terms of the GPLv3 license 
                       (see Amorpheus.py or LICENSE for details)
      """)

print('/\\'*32)
print()


if FileFormat == 3 or FileFormat == 4:   
    lambd= 1

# Names written in the parameters menu [Par]
ParListNames = [
'Filename',
'Filename background',
'File Format',
'incoherent (1 on empty cell, 0 on sample)',
'lambda (nm)',
'qmin (nm-1)',
'qmax (nm-1)',
'rmin (nm)',
'bkfactor',
'density (at/nm3)',
'Normalization Type',
'SofQ',
'QStep',
'Normalization save filename',
'Density rmin save filename',
'Density qmax save filename',
'Elements',
'Content']

#Real variables names
ParListRef = [
'filename',
'filename_bkg',
'FileFormat',
'incoherent',
'lambd',
'qmin',
'qmax',
'rmin',
'bkfactor',
'rhomi',
'Normalization',
'SofQ',
'QStep',
'F3_save_filename',
'F4_save_filename',
'F5_save_filename',
'elements',
'ratio']

ParListType = [
str,
str,
int,
int,
float,
float,
float,
float,
float,
float,
int,
int,
int,
str,
str,
str,
list,
np.ndarray]

FileFormat_names = [\
                    'not used',\
                    'Th No Background',\
                    'Th + Background ',\
                    'Q No Background',\
                    'Q + Background']
if Input_File:
    try:
        import configparser
        import ast
        
        def load_init_file(filename,ParListNames,ParListRef, ParListType):
            config = configparser.ConfigParser()
            config.read(filename)
            Names = ParListNames[:-2]
            Type = ParListType[:-2]
            
            var = [config['Parameters'][name] for name in Names]
            var_type = [t(x) for t,x in zip(Type,var)]
            for i,name in enumerate(ParListNames[:-2]):
                try:
                    exec(f'globals()["{ParListRef[i]}"] = {var_type[i]}')
                except:
                    exec(f'globals()["{ParListRef[i]}"] = "{var_type[i]}"')
            #special case of list of string
            var_type = [n.strip() for n in ast.literal_eval(config['Parameters'][ParListNames[-2]])]
            exec(f'globals()["{ParListRef[-2]}"] = {var_type}')
            var = config['Parameters'][ParListNames[-1]]
            var_type = np.array(var.replace('[','').replace(']','').split()).astype(np.float)
            exec(f'globals()["{ParListRef[-1]}"] = var_type')
        print()
        print(f" Settings from '{Init_filename}' file " + time.strftime("%Y/%m/%d %H:%M:%S"))
        load_init_file(Init_filename,ParListNames,ParListRef, ParListType)
        print()
    except:
        print()
        print(f"ERROR: No input file named '{Init_filename}' or corrupted file, using default settings from the script")
        print()
        print('To create a new Init file go in Parameters [p] and then [i] to create Init.txt')
        print()

else:
    print()
    print()
    print("Settings from the script, don't mess with the code !!!!")
    print()
    allPar = [eval(i) for i in ParListRef]
    eq.printPar(ParListNames,allPar)
    

elements,ratio = eq.ratio100(elements,ratio) #check sum = 1


F3_save_filename_ref = F3_save_filename
F4_save_filename_ref = F4_save_filename
F5_save_filename_ref = F5_save_filename

F3_save_pointer = int(0)
F4_save_pointer = int(0)
F5_save_pointer = int(0)


initialinitialrhomi = rhomi

thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi


while True:
    try:
        print('################################################################')
        print('#')      
        print('#  Visualize [v] - Normalize [n] - Loop rmin [r] - Loop qmax [x]\n#  Batch [b] - Chi2 [c] - Parameters [p] - About [a] -  Quit [q]')
        print('#')  
        print('################################################################')
        print()
        case = input()
        print()
        if case in ('v'):
            try:
                qmin, qmax, qmin_file, qmax_file = eq.plotFile(filename, filename_bkg, FileFormat, lambd, qmin, qmax)
                print('\nVisualization selected')
                print()
            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print() 
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)

    ############################################################
    
        elif case in ('n'):
            try:

                rmax=1
                
                print('\nNormalization selected')
                print()
                print(f'Your file is artificially cut from {qmin:.3f} to {qmax:.3f} nm^{-1}')
                print()
                print(f'average density: {rhomi} atoms/nm^3')
                print()
                [print(f'{el} content: {ratio[index]},', end = ' ') for index, el in enumerate(elements)]
                print()
                print(f'\nbkfactor: {bkfactor}, lambda: {lambd} nm')
                

                
                if FileFormat == 2:
                    th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio) #diff is just needed to define the values of S at low q values
                    Imeas = np.array(Imeas)
                    Ibkg = np.array(Ibkg)
                    Ibrut=Imeas-bkfactor*Ibkg
                
                elif FileFormat == 4:
                    q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio) 
                    Imeas = np.array(Imeas)
                    Ibkg = np.array(Ibkg)
                    Ibrut=Imeas-bkfactor*Ibkg
                    
                elif FileFormat == 3:
                    q,Ibrut, incohP, f2med, fmed2=eq.openFiles_3(filename,qmin,qmax,incoherent,elements,ratio) 
                    
                elif FileFormat == 1:
                    th,q,Ibrut, incohP, f2med, fmed2=eq.openFiles_1(filename,lambd,thmin,thmax,QStep,incoherent,elements,ratio) 
        
        
        
                if FileFormat == 2 or FileFormat == 4:
                    fig_title =  'Intensity'
                    plt.figure(fig_title)
                    plt.plot(q,Imeas,label='I$^{meas}$')
                    plt.plot(q,Ibkg,label='I$^{bkg}$')
                    plt.plot(q,Ibrut,label='I$^{sample}$')
                    plt.xlabel('q (nm$^{-1}$)')
                    plt.ylabel('I (a.u.)')
                    plt.grid()
                    plt.legend(loc='best')
                    plt.title(fig_title)
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                    
                    
                elif FileFormat == 3 or FileFormat == 1:
                    fig_title =  'Intensity'
                    plt.figure(fig_title)
                    plt.plot(q,Ibrut,label='I$^{sample}$')
                    plt.xlabel('q (nm$^{-1}$)')
                    plt.ylabel('I (a.u.)')
                    plt.legend(loc='best')
                    plt.grid()
                    plt.legend(loc='best')
                    plt.title(fig_title)
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                              
                    
                R=np.arange(0,rmax,0.0005)
                R[0]=0.000001 
                        
                R=np.linspace(0,rmax,9001)
                R[0]=5.00000000e-04
                q_old=q
                pasq = q_old[1]-q_old[0]   
                addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                q=np.concatenate((addq,q_old))
                
                S = np.zeros((eq.prout+1,len(q)))
                gr = np.zeros((eq.prout+1,len(R)))
                F = np.zeros((eq.prout+1,len(R)))
                delF = np.zeros((len(R)))
                Sinv = np.zeros((len(q)))
                
                diff = np.abs(q-q_old[0]).argmin()
                S,F,gder = eq.getSandF3(q,S,F,gr,Ibrut,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2, SofQ, delF, Sinv)
                
                fig_title =  'Structure factor S(q)'
                plt.figure(fig_title)
                plt.title(fig_title)
                [ plt.plot(q,S[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                plt.legend()
                plt.xlabel('q(nm$^{-1}$)')
                plt.ylabel('S(q)')
                plt.grid()        
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode
                
                fig_title =  'Distribution function F(r)'
                plt.figure(fig_title)
                plt.title(fig_title)
                [ plt.plot(R,F[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                plt.legend()
                plt.xlabel('r(nm)')
                plt.ylabel('F(r)')
                plt.grid()          
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode
                   
                
                indrmin = (np.abs(R - rmin)).argmin()
                dens=np.zeros(eq.prout+1)
                for i in range(0,eq.prout+1):
                    slope, intercept, r_value, p_value, std_err = stats.linregress(R[0:indrmin], F[i][0:indrmin])
                    dens[i] = -slope/(4*np.pi)
                
                if debug:
                    fig_title =  'Density as linear slope'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    plt.scatter(np.arange(eq.prout+1),dens)
                    plt.xlabel('prout')
                    plt.ylabel('density as lin slope (at/nm3)')
                    plt.grid()          
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                
                print()
    
                if debug:    
                    print(f"slope: {slope:.3f}, intercept: {intercept:.3f}, std_err: {std_err:.3f}")  
                    print()  
                    print(f'The density from the linear fit of F(r) at low r: {dens[eq.prout-1]}')
                    print()  
                
                
                intmin=rmin
                intMin=np.abs(R-intmin).argmin()
        
        
                for index in range(0,len(gder)):
                    peaks, _  = find_peaks(gder[index])
                    intMax = peaks[np.argmax(peaks>intMin)]
                    CN1=2*np.trapz(4*np.pi*rhomi*R[intMin:intMax]**2*gder[index][intMin:intMax],R[intMin:intMax])
                    peak_2 = peaks[np.argmax(peaks>intMax+1)]
                    intMax2=intMax+gder[index][intMax:peak_2].argmin()
                    CN2=np.trapz(4*np.pi*rhomi*R[intMin:intMax2]**2*gder[index][intMin:intMax2],R[intMin:intMax2])
                    if debug:
                        fig_title = f'debug_{index}'
                        plt.figure(fig_title)
                        plt.plot(R,gder[index])
                        plt.plot(R[peaks],gder[index][peaks],'or')
                        plt.plot(R[intMax2],gder[index][intMax2],'ob')
                        plt.axvline(R[intMin],label = f'rmin_i = {intMin}, val = {R[intMin]:.3f}')
                        plt.axvline(R[intMax],label = f'CN1_i = {intMax}, val = {R[intMax]:.3f}')
                        plt.axvline(R[intMax2],label = f'CN2_i = {intMax2}, val = {R[intMax2]:.3f}')
                        plt.title(f'g(r) {index}')
                        plt.legend()
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                    print(f'{index}: CN1 = {CN1:.2f}, CN2 = {CN2:.2f}\tr0={R[intMin]:.3f}, rmax={R[intMax]:.3f}, rint={R[intMax2]:.3f}')
                
                fig_title =  'Radial distribution function g(r)'
                plt.figure(fig_title)
                plt.title(fig_title)
                [ plt.plot(R,gder[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                plt.legend()
                plt.xlabel('r(nm)')
                plt.ylabel('g(r)') 
                plt.grid()
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode
                
                if Normalization ==1:
                    happy = 0
                    while happy == 0:            
                        modfa = float(input(f"Choose parameter 'a' as the position of the first peak of the g(r) in nm ({R[intMax]:.3f}): "))
                        modfb = float(input(f"Please chose parameter 'b' as the width of the first peak of the g(r) in nm ({R[intMax]-rmin:.3f}): "))
                        
                        fig_title =  'Radial distribution function g(r)'
                        plt.figure(fig_title)
                        plt.title(fig_title)
                        plt.plot(R,gr[-1],label=('g(r)'))
                        Flast = eq.calcF_speed_ModificationFunction(R,q,S[-1],modfa,modfb)
                        grlast=1+Flast/(4*np.pi*R*rhomi)
                        plt.plot(R,grlast,label=('g(r) with M(Q,$\Delta$(r))'),color='k')
                        deltaR=np.pi/q[-1]*(1-np.exp(-np.abs(R[:,None]-modfa)/modfb))
                        normlzz=np.max(gr[-1])/np.max(deltaR)
                        plt.plot(R,normlzz*deltaR,'--',label=f'$\Delta$(r)*{int(normlzz)}: a={modfa}nm\nb={modfb}nm')
                        plt.legend()
                        plt.xlabel('r(nm)')
                        plt.ylabel('g(r)')
                        plt.legend(loc='best')
                        np.savetxt(F3_save_filename+'_grlast.txt',np.vstack((R,grlast)).T)
                        np.savetxt(F3_save_filename+'_deltaR.txt',np.vstack((R,deltaR.T)).T)
                        plt.savefig(F3_save_filename+'ModificationFunction.png',bbox_inches='tight',dpi=DPI,format='png')
                        plt.savefig(F3_save_filename+'ModificationFunction.eps',bbox_inches='tight',dpi=DPI,format='eps')                
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                        print(f'g(r) with a={modfa} nm b={modfb} nm')
                        
                        happy = float(input('If you like the result choose 1, otherwise 0: '))

                print()
                print('Save ? y/n')
                case = input()
                if case in ('Y','y'):
                    F3_save_filename, F3_save_filename_ref, F3_save_pointer = eq.rename_F3(F3_save_filename, F3_save_filename_ref, F3_save_pointer)
                    if FileFormat == 2 or FileFormat == 4:
                        np.savetxt(F3_save_filename+'_I.txt',np.vstack((q_old,Imeas,Ibkg,Ibrut)).T)
                    elif FileFormat == 3 or FileFormat == 1:
                        np.savetxt(F3_save_filename+'_I.txt',np.vstack((q_old,Ibrut)).T)
                    np.savetxt(F3_save_filename+'_S.txt',np.vstack((q,S)).T)
                    np.savetxt(F3_save_filename+'_S_lastcolumn.txt',np.vstack((q,S[(len(S)-1)])).T)
                    np.savetxt(F3_save_filename+'_F.txt',np.vstack((R,F)).T)
                    np.savetxt(F3_save_filename+'_gder.txt',np.vstack((R,gder)).T)
                    
                    fig_title =  'Intensity'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    if FileFormat == 2 or FileFormat == 4:
                        plt.plot(q_old,Imeas, label='I$^{meas}$')
                        plt.plot(q_old,Ibkg,label='I$^{bkg}$')
                    plt.plot(q_old,Ibrut, label='I$^{sample}$')
                    plt.xlabel('q (nm$^{-1}$)')
                    plt.ylabel('I (a.u.)')
                    plt.legend(loc='best')
                    plt.grid()

                    plt.savefig(F3_save_filename+'_I',bbox_inches='tight',dpi=DPI)
                    plt.close()
                    
                    fig_title =  'Structure factor S(q)'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    [ plt.plot(q,S[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                    plt.legend()
                    plt.xlabel('q(nm$^{-1}$)')
                    plt.ylabel('S(q)')  
                    plt.grid()
                    plt.savefig(F3_save_filename+'_S',bbox_inches='tight',dpi=DPI)
                    plt.close()
                    
                    fig_title =  'Distribution function F(r)'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    [ plt.plot(R,F[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                    plt.legend()
                    plt.xlabel('r(nm)')
                    plt.ylabel('F(r)')  
                    plt.grid()

                    plt.savefig(F3_save_filename+'_F',bbox_inches='tight',dpi=DPI)
                    plt.close()
                    
                    fig_title =  'Radial distribution function g(r)'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    [ plt.plot(R,gder[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                    plt.legend()
                    plt.xlabel('r(nm)')
                    plt.ylabel('g(r)')  
                    plt.grid()

                    plt.savefig(F3_save_filename+'_gder',bbox_inches='tight',dpi=DPI)
                    plt.close()

                    if Normalization ==1:        
                        fig_title =  'Radial distribution function g(r)'
                        plt.figure(fig_title)
                        plt.title(fig_title)
                        plt.plot(R,gr[-1],label=('g(r)'))
                        Flast = eq.calcF_speed_ModificationFunction(R,q,S[-1],modfa,modfb)
                        grlast=1+Flast/(4*np.pi*R*rhomi)
                        plt.plot(R,grlast,label=('g(r) with M(Q,$\Delta$(r))'),color='k')
                        deltaR=np.pi/q[-1]*(1-np.exp(-np.abs(R[:,None]-modfa)/modfb))
                        normlzz=np.max(gr[-1])/np.max(deltaR)
                        plt.plot(R,normlzz*deltaR,'--',label=f'$\Delta$(r)*{int(normlzz)}: a={modfa}nm\nb={modfb}nm')
                        plt.xlabel('r(nm)')
                        plt.ylabel('g(r)')
                        plt.legend(loc='best')
                        np.savetxt(F3_save_filename+'_grlast.txt',np.vstack((R,grlast)).T)
                        np.savetxt(F3_save_filename+'_deltaR.txt',np.vstack((R,deltaR.T)).T)
                        plt.savefig(F3_save_filename+'ModificationFunction.png',bbox_inches='tight',dpi=DPI,format='png')
                        plt.savefig(F3_save_filename+'ModificationFunction.eps',bbox_inches='tight',dpi=DPI,format='eps')
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode

            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print() 
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)
            
    ############################################################
    
    
            
        elif case in ('F4rough','f4rough','r', 'R'):
            try:


                rmax=0.4 
                
                print('\nDensity rmin selected')
                print()
                print(f'average density: {rhomi} atoms/nm^3')
                [print(f'{el} content: {ratio[index]},', end = ' ') for index, el in enumerate(elements)] 
                print(f'\nbkfactor: {bkfactor}, lambda: {lambd} nm')
                             
                    
                if FileFormat == 2:
                    th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio) #diff is just needed to define the values of S at low q values
                    Imeas = np.array(Imeas)
                    Ibkg = np.array(Ibkg)
                    Ibrut=Imeas-bkfactor*Ibkg
                
                elif FileFormat == 4:
                    q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio) 
                    Imeas = np.array(Imeas)
                    Ibkg = np.array(Ibkg)
                    Ibrut=Imeas-bkfactor*Ibkg
                    
                elif FileFormat == 3:
                    q,Ibrut, incohP, f2med, fmed2=eq.openFiles_3(filename,qmin,qmax,incoherent,elements,ratio) 
                    
                elif FileFormat == 1:
                    th,q,Ibrut, incohP, f2med, fmed2=eq.openFiles_1(filename,lambd,thmin,thmax,QStep,incoherent,elements,ratio)
        
        
                
                
                print('\nEnter init rhomi (at/nm^3)')
                try: initrhomi = float(input())
                except Exception as e: print(e)
                print('\ninit rhomi = '+ str(initrhomi))
         
                if FileFormat == 2 or FileFormat == 4:
                    print('\nEnter init bkfactor')
                    try: initbkfactor = float(input())
                    except Exception as e: print(e)
                    print('\ninit bkfactor = '+ str(initbkfactor))
                
                fit_params = Parameters()
                fit_params.add('rhomi', value = initrhomi, min = max(1,initrhomi-50), max = initrhomi+50)
                if FileFormat == 2 or FileFormat == 4:
                    fit_params.add('bkfactor', value = initbkfactor, min = 0.7, max = 1.2)
                
                
                
                R=np.arange(0,rmax,0.0005)
                R[0]=0.000001 
                q_old=q
                pasq = q_old[1]-q_old[0]   
                addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                q=np.concatenate((addq,q_old))
                
                S = np.zeros((eq.prout+1,len(q)))
                gr = np.zeros((eq.prout+1,len(R)))
                F = np.zeros((eq.prout+1,len(R)))
                delF = np.zeros((len(R)))
                Sinv = np.zeros((len(q)))

                rmin_min, rmin_max, steprmin, array_rmin = eq.input_rmin()
                
                results = np.zeros((len(array_rmin),4))
                diff = np.abs(q-q_old[0]).argmin()
                
                F4_save_filename, F4_save_filename_ref, F4_save_pointer = eq.rename_F45(F4_save_filename, F4_save_filename_ref, F4_save_pointer, filename, filename_bkg, rmin_min, rmin_max)
                
                print()
                print('Control+C to interrupt')
                print()
                if FileFormat == 2 or FileFormat == 4:
                    interrupt = 0
                    try:
                        for i, _rmin in enumerate(array_rmin):
                            print(f'rmin = {_rmin}, loop {i} / {len(array_rmin)-1}' )
                            ind = (np.abs(R - _rmin)).argmin()
                            now = time.time()
                            args = (q,S,F,gr,Imeas,Ibkg,R,_rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                            out = minimize(eq.chisquare_DAC_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})
                            out_rhomi = out.params["rhomi"].value
                            out_bkfactor = out.params["bkfactor"].value
                            chi2function = out.residual[0]
                            print(f'rmin = {_rmin:.3f} rhomi = {out_rhomi:.3f} bkfactor = {out_bkfactor:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                            results[i] = [_rmin, out_rhomi, out_bkfactor, chi2function]
                            interrupt = interrupt +1
                    except KeyboardInterrupt:
                        results = results[0:interrupt,:]
                        print()
                        print('Interrupted !')
                        print()
                                
                    
                elif FileFormat == 3 or FileFormat == 1:
                    interrupt = 0
                    try:              
                        for i, _rmin in enumerate(array_rmin):
                            print(f'rmin = {_rmin}, loop {i} / {len(array_rmin)-1}' )
                            ind = (np.abs(R - _rmin)).argmin()
                            now = time.time()
                            args = (q,S,F,gr,Ibrut,R,_rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                            out = minimize(eq.chisquare_LV_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})
                            out_rhomi = out.params["rhomi"].value
                            chi2function = out.residual[0]
                            print(f'rmin = {_rmin:.3f} rhomi = {out_rhomi:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                            results[i] = [_rmin, out_rhomi, 0,chi2function]
                                
                            interrupt = interrupt +1
                    except KeyboardInterrupt:
                        results = results[0:interrupt,:]
                        print()
                        print('Interrupted !')
                        print()
                        
                allPar = [eval(i) for i in ParListRef]
                header = eq.header_r(allPar , rmin_min, rmin_max, initialinitialrhomi,FileFormat_names[FileFormat])
                
                """ATTENTION: if change name in savetxt change name in eq.rename_F45 too !!!!"""
                
                np.savetxt(f'{F4_save_filename}_{filename[0:-4]}-{filename_bkg[0:-3]}_{rmin_min}-{rmin_max}.txt',results, header = header)
                print(f'The file with the results has been saved as: {F4_save_filename}_{filename[0:-4]}-{filename_bkg}_{rmin_min}-{rmin_max}.txt')
        
                
                if FileFormat == 2 or FileFormat == 4:
                    
                    fig_title =  'Loop over rmin'
                    fig, ax1 = plt.subplots()
                    ax1.set_title(fig_title)
                    color = 'tab:red'
                    ax1.set_xlabel('rmin (nm)')
                    ax1.set_ylabel('bkfactor', color=color)
                    ax1.plot(np.transpose(results)[0], np.transpose(results)[2], color=color)
                    ax1.tick_params(axis='y', labelcolor=color)
                    ax2 = ax1.twinx()
                    color = 'tab:blue'
                    ax2.set_ylabel('chisquare', color=color)
                    ax2.plot(np.transpose(results)[0], np.transpose(results)[3], color=color)
                    ax2.tick_params(axis='y', labelcolor=color)
                    fig.savefig(F4_save_filename+'bkfactor'+str(rmin_min)+'_'+str(rmin_max)+'.png',dpi=DPI)       
                    plt.title(fig_title)
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                
                fig_title =  'Loop over rmin'
                fig, ax1 = plt.subplots()
                ax1.set_title(fig_title)
                color = 'tab:red'
                ax1.set_xlabel('rmin (nm)')
                ax1.set_ylabel('rho (at/nm$^3$)', color=color)
                ax1.plot(np.transpose(results)[0], np.transpose(results)[1], color=color)
                ax1.tick_params(axis='y', labelcolor=color)
                ax2 = ax1.twinx()
                color = 'tab:blue'
                ax2.set_ylabel('chisquare', color=color)
                ax2.plot(np.transpose(results)[0], np.transpose(results)[3], color=color)
                ax2.tick_params(axis='y', labelcolor=color)
                fig.savefig(F4_save_filename+'_chisquare_'+str(rmin_min)+'_'+str(rmin_max)+'.png',dpi=DPI)
                plt.title(fig_title)
                plt.show() 
                plt.pause(0.001) #used to update figure in terminal mode
                
            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print() 
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)
            
    ############################################################
            
        elif case in ('Batch','batch','b', 'B'):
            try:

                rmax=1
                
                print('\nBatch selected')
                print()
                batch_data, temperatures, pressures, folders, batch_loop_n = eq.batch_check()
                print(f"{batch_loop_n} files from Batch.ods")
                print()
                
                timetosave = time.strftime("%Y%m%d_%Hh%M")
                foldertosavebatch = f'Batch_{timetosave}'
                if not os.path.exists(foldertosavebatch):
                    os.makedirs(foldertosavebatch)
                      
                
                Names = ParListNames[:-5]
                Type = ParListType[:-5]

                
                results = np.zeros((batch_loop_n,4))
                
                S_batch = []
                q_batch = []
                gder_batch = []
                r_batch = []
                glast_batch = []
                rlast_batch = []

                try:
                    interrupt = 0
                    for batch_line in range(batch_loop_n):
                        
                        
                        var = [i for i in batch_data[1+batch_line,3:16] if i != None]
                        var_type = [t(x) for t,x in zip(Type,var)]
                        
                        for i,name in enumerate(ParListNames[:-5]):
                            try:
                                exec(f'globals()["{ParListRef[i]}"] = {var_type[i]}')
                            except:
                                exec(f'globals()["{ParListRef[i]}"] = "{var_type[i]}"')
                        
                        #special case of list of string
                        var_type = [batch_data[1+batch_line,16],batch_data[1+batch_line,18],batch_data[1+batch_line,20],batch_data[1+batch_line,22],batch_data[1+batch_line,24],batch_data[1+batch_line,26]]
                        exec(f'globals()["{ParListRef[-2]}"] = {var_type}')

                        
                        var_type = np.array([batch_data[1+batch_line,16+1],batch_data[1+batch_line,18+1],batch_data[1+batch_line,20+1],batch_data[1+batch_line,22+1],batch_data[1+batch_line,24+1],batch_data[1+batch_line,26+1]])
                        exec(f'globals()["{ParListRef[-1]}"] = var_type')

                            
                        print()
                        print(f'Batch loop {batch_line+1} over {batch_loop_n}')
                        print()
                        
                        inputfilename = filename
                        if filename[-2:] == 'xt' or filename[-2:] == 'at':
                            inputfilename = filename[:-4]
                        elif filename[-2:] == 'xy':
                            inputfilename = filename[:-3]
                            
                        
                        try:
                            if (Path(folders[batch_line]) / filename).is_file():
                                #if no subfolder and files in the same folder
                                filename = folders[batch_line] / filename
                                filename_bkg = folders[batch_line] / filename_bkg
                            else:
                                print('Batch file not found')
                        except Exception:
                                #mix of space and empty
                            pass
            
                        
                            
                        print('Batch loop parameters:')
                        print()
                        allPar = [eval(i) for i in ParListRef]
                        eq.printPar(ParListNames,allPar)
                        print()
                                     
                            
                        if FileFormat == 2:
                            thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                            thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi
                            th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio) 
                            Imeas = np.array(Imeas)
                            Ibkg = np.array(Ibkg)
                            Ibrut=Imeas-bkfactor*Ibkg
                        
                        elif FileFormat == 4:
                            q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio) 
                            Imeas = np.array(Imeas)
                            Ibkg = np.array(Ibkg)
                            Ibrut=Imeas-bkfactor*Ibkg
                            
                        elif FileFormat == 3:
                            q,Ibrut, incohP, f2med, fmed2=eq.openFiles_3(filename,qmin,qmax,incoherent,elements,ratio) 
                            
                        elif FileFormat == 1:
                            thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                            thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi
                            th,q,Ibrut, incohP, f2med, fmed2=eq.openFiles_1(filename,lambd,thmin,thmax,QStep,incoherent,elements,ratio)
                
                
                        
                        fit_params = Parameters()
                        fit_params.add('rhomi', value = rhomi, min = max(1,rhomi-50), max = rhomi+50)
                        if FileFormat == 2 or FileFormat == 4:
                            fit_params.add('bkfactor', value = bkfactor, min = 0.7, max = 1.2)
                        
                        
                        
                        R=np.arange(0,rmax,0.0005)
                        R[0]=0.000001
                        q_old=q
                        pasq = q_old[1]-q_old[0]   
                        addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                        q=np.concatenate((addq,q_old))
                        
                        S = np.zeros((eq.prout+1,len(q)))
                        gr = np.zeros((eq.prout+1,len(R)))
                        F = np.zeros((eq.prout+1,len(R)))
                        delF = np.zeros((len(R)))
                        Sinv = np.zeros((len(q)))

                        

                        diff = np.abs(q-q_old[0]).argmin()
                        
                        print('Fitting ...')
                        print()
                        
                        if FileFormat == 2 or FileFormat == 4:
                            
                            ind = (np.abs(R - rmin)).argmin()
                            now = time.time()
                            args = (q,S,F,gr,Imeas,Ibkg,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                            out = minimize(eq.chisquare_DAC_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})
                            out_rhomi = out.params["rhomi"].value
                            out_bkfactor = out.params["bkfactor"].value
                            chi2function = out.residual[0]
                            print(f'rmin = {rmin:.3f} rhomi = {out_rhomi:.3f} bkfactor = {out_bkfactor:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                            results[batch_line] = [rmin, out_rhomi, out_bkfactor, chi2function]
                            rhomi = out_rhomi
                        

                        elif FileFormat == 3 or FileFormat == 1:
          
                            ind = (np.abs(R - rmin)).argmin()
                            now = time.time()
                            args = (q,S,F,gr,Ibrut,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                            out = minimize(eq.chisquare_LV_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})
                            out_rhomi = out.params["rhomi"].value
                            chi2function = out.residual[0]
                            print(f'rmin = {rmin:.3f} rhomi = {out_rhomi:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                            results[batch_line] = [rmin, out_rhomi, 0,chi2function]
                            rhomi = out_rhomi
                                

                        
                        if FileFormat == 3 or FileFormat == 1:
                            bkfactor = 0
                        else:
                            bkfactor = out_bkfactor
                        
            ###########c & p [x]
                                
                        print('\nNormalization')
                        print()
                        print(f'rhomi = {out_rhomi:.3f}')
            
                        if FileFormat == 2:

                            thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                            thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi
                            th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio) 
                            Imeas = np.array(Imeas)
                            Ibkg = np.array(Ibkg)
                            Ibrut=Imeas-bkfactor*Ibkg
                        
                        elif FileFormat == 4:
                            q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio) 
                            Imeas = np.array(Imeas)
                            Ibkg = np.array(Ibkg)
                            Ibrut=Imeas-bkfactor*Ibkg
                            
                        elif FileFormat == 3:
                            q,Ibrut, incohP, f2med, fmed2=eq.openFiles_3(filename,qmin,qmax,incoherent,elements,ratio) 
                            
                        elif FileFormat == 1:
                            thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                            thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi
                            th,q,Ibrut, incohP, f2med, fmed2=eq.openFiles_1(filename,lambd,thmin,thmax,QStep,incoherent,elements,ratio)
                        
                        fig_title =  f'I {inputfilename} {pressures[batch_line]} GPa'

                        
                        if FileFormat == 2 or FileFormat == 4:
                            
                            plt.figure(fig_title)
                            plt.plot(q,Imeas,label='I$^{meas}$')
                            plt.plot(q,Ibkg,label='I$^{bkg}$')
                            plt.plot(q,Ibrut,label='I$^{sample}$')
                            plt.xlabel('q (nm$^{-1}$)')
                            plt.ylabel('I (a.u.)')
                            plt.grid()
                            plt.legend(loc='best')
                            plt.title(fig_title)
                            plt.savefig(f'{foldertosavebatch}/I_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                            plt.show()
                            plt.pause(0.001) #used to update figure in terminal mode
                            
                            
                        elif FileFormat == 3 or FileFormat == 1:
                            plt.figure(fig_title)
                            plt.plot(q,Ibrut,label='I$^{sample}$')
                            plt.xlabel('q (nm$^{-1}$)')
                            plt.ylabel('I (a.u.)')
                            plt.legend(loc='best')
                            plt.grid()
                            plt.legend(loc='best')
                            plt.title(fig_title)
                            plt.savefig(f'{foldertosavebatch}/I_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                            plt.show()
                            plt.pause(0.001) #used to update figure in terminal mode
                                  
                            
                
                        R=np.linspace(0,rmax,9001)
                        R[0]=5.00000000e-04
                        q_old=q
                        pasq = q_old[1]-q_old[0]   
                        addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                        q=np.concatenate((addq,q_old))
                        
                        S = np.zeros((eq.prout+1,len(q)))
                        gr = np.zeros((eq.prout+1,len(R)))
                        F = np.zeros((eq.prout+1,len(R)))
                        delF = np.zeros((len(R)))
                        Sinv = np.zeros((len(q)))
                        
                        diff = np.abs(q-q_old[0]).argmin()

                        S,F,gder = eq.getSandF3(q,S,F,gr,Ibrut,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2, SofQ, delF, Sinv)
                        
                        fig_title = f'S(q) {inputfilename} {pressures[batch_line]} GPa'

                        plt.figure(fig_title)
                        plt.title(fig_title)
                        [ plt.plot(q,S[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                        plt.legend()
                        plt.xlabel('q(nm$^{-1}$)')
                        plt.ylabel('S(q)')
                        plt.grid()        

                        plt.title(fig_title)
                        plt.savefig(f'{foldertosavebatch}/S_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                        tosave=np.vstack([q,S])
                        np.savetxt(f'{foldertosavebatch}/S_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.txt',tosave.T)
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
 
                        
                        S_batch.append(S[eq.prout])
                        q_batch.append(q)
                        
                        fig_title = f'F(r) {inputfilename} {pressures[batch_line]} GPa'

                        plt.figure(fig_title)
                        [ plt.plot(R,F[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                        plt.legend()
                        plt.xlabel('r(nm)')
                        plt.ylabel('F(r)')
                        plt.grid()
                        plt.title(fig_title)
                        plt.savefig(f'{foldertosavebatch}/F_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                        tosave=np.vstack([R,F])
                        np.savetxt(f'{foldertosavebatch}/F_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.txt',tosave.T)                        
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                        
                                
                        
                        indrmin = (np.abs(R - rmin)).argmin()
                        dens=np.zeros(eq.prout+1)
                        for i in range(0,eq.prout+1):
                            slope, intercept, r_value, p_value, std_err = stats.linregress(R[0:indrmin], F[i][0:indrmin])
                            dens[i] = -slope/(4*np.pi)

                        print(f'The density from linear fit of the F(r) for r < rmin is: {dens[eq.prout-1]:.3f}')
                        print()

                        intmin=rmin
                        intMin=np.abs(R-intmin).argmin()
            
                        for index in range(0,len(gder)):
                            peaks, _  = find_peaks(gder[index])
                            intMax = peaks[np.argmax(peaks>intMin)]
                            CN1=2*np.trapz(4*np.pi*rhomi*R[intMin:intMax]**2*gder[index][intMin:intMax],R[intMin:intMax])
                            peak_2 = peaks[np.argmax(peaks>intMax+1)]
                            intMax2=intMax+gder[index][intMax:peak_2].argmin()
                            CN2=np.trapz(4*np.pi*rhomi*R[intMin:intMax2]**2*gder[index][intMin:intMax2],R[intMin:intMax2])
                            if debug:
                                fig_title = f'debug_{index}'
                                plt.figure(fig_title)
                                plt.plot(R,gder[index])
                                plt.plot(R[peaks],gder[index][peaks],'or')
                                plt.plot(R[intMax2],gder[index][intMax2],'ob')
                                plt.axvline(R[intMin],label = f'rmin_i = {intMin}, val = {R[intMin]:.3f}')
                                plt.axvline(R[intMax],label = f'CN1_i = {intMax}, val = {R[intMax]:.3f}')
                                plt.axvline(R[intMax2],label = f'CN2_i = {intMax2}, val = {R[intMax2]:.3f}')
                                plt.title(f'g(r) {index}')
                                plt.legend()
                                plt.show()
                                plt.pause(0.001) #used to update figure in terminal mode
                            print(f'{index}: CN1 = {CN1:.2f}, CN2 = {CN2:.2f}\tr0={R[intMin]:.3f}, rmax={R[intMax]:.3f}, rint={R[intMax2]:.3f}')
                        
                        fig_title = f'g(r) {inputfilename} {pressures[batch_line]} GPa'

                        plt.figure(fig_title)
                        plt.title(fig_title)
                        [ plt.plot(R,gder[i],label=(f'{i}')) for i in range(0,eq.prout+1)]
                        plt.legend()
                        plt.xlabel('r(nm)')
                        plt.ylabel('g(r)') 
                        plt.grid()
                        plt.savefig(f'{foldertosavebatch}/gder_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                        tosave=np.vstack([R,gder])
                        np.savetxt(f'{foldertosavebatch}/gder_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.txt',tosave.T)                        
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                        
                        gder_batch.append(gder[eq.prout])
                        r_batch.append(R)
                        
                        modfa = R[intMax]
                        modfb = modfa - rmin
                        
                        if Normalization ==1:
                            fig_title = f'g(r) {inputfilename} {pressures[batch_line]} GPa'

                            plt.figure(fig_title)
                            plt.plot(R,gr[-1],label=('g(r)'))
                            Flast = eq.calcF_speed_ModificationFunction(R,q,S[-1],modfa,modfb)
                            grlast=1+Flast/(4*np.pi*R*rhomi)
                            plt.plot(R,grlast,label=('g(r) with M(Q,$\Delta$(r))'),color='k')
                            deltaR=np.pi/q[-1]*(1-np.exp(-np.abs(R[:,None]-modfa)/modfb))
                            normlzz=np.max(gr[-1])/np.max(deltaR)
                            plt.plot(R,normlzz*deltaR,'--',label=f'$\Delta$(r)*{normlzz:.2f}: a={modfa:.3f}nm\nb={modfb:.3f}nm')
                            plt.legend()
                            plt.xlabel('r(nm)')
                            plt.ylabel('g(r)')

                            plt.savefig(f'{foldertosavebatch}/grlast_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.png',dpi=DPI)
                            tosave=np.array([R,grlast])
                            np.savetxt(f'{foldertosavebatch}/grlast_{inputfilename}_{pressures[batch_line]}'.replace('.','p')+'.txt',tosave.T)                            
                            plt.show()
                            plt.pause(0.001) #used to update figure in terminal mode
                            print(f'g(r) with a={modfa:.3f} nm b={modfb:.3f} nm')
                            glast_batch.append(grlast)
                            rlast_batch.append(R)
                                    
                                                                
                        interrupt = interrupt +1 
                except KeyboardInterrupt:
                   results = results[0:interrupt,:]
                   pressures = pressures[0:interrupt]
                   temperatures = temperatures[0:interrupt]
                   print()
                   print('Interrupted !')
                   print()   
                trasl=0.6
                
                fig_title =  'Batch results S(q)'

                plt.figure(fig_title)
                plt.title(fig_title)
                for i in range(len(S_batch)):
                    lab = (f'{pressures[i]} GPa, {temperatures[i]} K')
                    plt.plot(q_batch[i],S_batch[i]+i*trasl, label = lab)
                plt.xlabel('q(nm$^{-1}$)')
                plt.ylabel('S(q)')
                plt.legend()
                plt.savefig(f'{foldertosavebatch}/Batch_results_S.png',dpi=DPI)
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode
                
                fig_title =  'Batch results g(r)'

                plt.figure(fig_title)
                plt.title(fig_title)
                for i in range(len(gder_batch)):
                    lab = (f'{pressures[i]} GPa, {temperatures[i]} K')
                    plt.plot(r_batch[i],gder_batch[i]+i*trasl, label = lab)
                plt.xlabel('r(nm)')
                plt.ylabel('g(r)')
                plt.legend()
                plt.savefig(f'{foldertosavebatch}/Batch_results_gder.png',dpi=DPI)
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode                
                
                if Normalization == 1:
                    fig_title =  'Batch results g(r)'
                    plt.figure(fig_title)
                    plt.title(fig_title)
                    for i in range(len(glast_batch)):
                        lab = (f'{pressures[i]} GPa, {temperatures[i]} K')
                        plt.plot(rlast_batch[i],glast_batch[i]+i*trasl, label = lab)
                    plt.xlabel('r(nm)')
                    plt.ylabel('g(r) with Modification Function')
                    plt.legend()
                    plt.savefig(f'{foldertosavebatch}/Batch_results_glast.png',dpi=DPI)
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                
                fig_title =  'Batch results g(r)'
                plt.figure(fig_title)
                plt.title(fig_title)
                plt.scatter(pressures,results[:,1])   
                plt.xlabel('P(GPa)')
                plt.ylabel('density (at/nm$^{3}$)')
                plt.savefig(f'{foldertosavebatch}/Batch_results_density.png',dpi=DPI)
                plt.show()
                plt.pause(0.001) #used to update figure in terminal mode
                
                headbatch = 'pressure, temperature, rmin, density, bkfactor, chisquare '
                riassunto=np.array([pressures,temperatures,results[:,0],results[:,1],results[:,2],results[:,3]])
                np.savetxt(f'{foldertosavebatch}/Batch_results_{timetosave}.txt',riassunto.T, header = headbatch, comments='')

                eq.batch_save(f'{foldertosavebatch}/Batch_{timetosave}.ods')

            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print()
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)
    ############################################################        
            
            
        elif case in ('X','x'):
            try:


                rmax=0.4
                
                print('\Density qmax selected')
                
                qmin, qmax, qmin_file, qmax_file = eq.checkFile(filename, filename_bkg, FileFormat, lambd, qmin, qmax)                

                print(f'average density: {rhomi} atoms/nm^3')
                print()
                [print(f'{el} content: {ratio[index]},', end = ' ') for index, el in enumerate(elements)] 
                print()
                print(f'\nbkfactor: {bkfactor}, lambda: {lambd} nm')   
                    
                
                print('\nEnter init rhomi (at/nm^3)')
                try: initrhomi = float(input())
                except Exception as e: print(e)
                print('\ninit rhomi = '+ str(initrhomi))
         
                if FileFormat == 2 or FileFormat == 4:
                    print('\nEnter init bkfactor')
                    try: initbkfactor = float(input())
                    except Exception as e: print(e)
                    print('\ninit bkfactor = '+ str(initbkfactor))
                
                fit_params = Parameters()
                fit_params.add('rhomi', value = initrhomi, min = max(1,initrhomi-50), max = initrhomi+50)
                if FileFormat == 2 or FileFormat == 4:
                    fit_params.add('bkfactor', value = initbkfactor, min = 0.7, max = 1.2)
        
        
        
                
                print(f'\nEnter qmax_min (nm^-1), min value larger than {qmin:.3f}')
                try: qmax_min = float(input())
                except Exception as e: print(e)
                print('\ninit qmax_min = '+ str(qmax_min))
                if qmax_min <= qmin : raise ValueError(('/\\'*30+"\n\nError, qmax_min must be larger than qmin, check Parameters [p]\n\n"+'/\\'*30))
                
                print(f'\nEnter qmax_max (nm^-1), max value lower than {qmax_file:.3f}')
                try: qmax_max = float(input())
                except Exception as e: print(e)
                print('\ninit qmax_max = '+ str(qmax_max))        
                if qmax_max >= qmax_file :
                    qmax_max = qmax_file
                    print(f'\nThe value you chose is not acceptable, qmax_max set to {qmax_file} nm^-1')
                
                print('\nEnter the step (nm^-1), suggested value 5 or 10')
                try: stepqmax = float(input())
                except Exception as e: print(e)

                array_qmax=eq.arange(qmax_min,qmax_max,stepqmax)
                
        
                
                F5_save_filename, F5_save_filename_ref, F5_save_pointer = eq.rename_F45(F5_save_filename, F5_save_filename_ref, F5_save_pointer, filename, filename_bkg, qmax_min, qmax_max)
                results = np.zeros((len(array_qmax),4))
                
                fit_params = Parameters()
                fit_params.add('rhomi', value = initrhomi, min = max(1,initrhomi-50), max = initrhomi+50)
                if FileFormat == 2 or FileFormat == 4:
                    fit_params.add('bkfactor', value = initbkfactor, min = 0.7, max = 1.2)
                    
        
                print()
                print('Control+C to interrupt')
                print()
                try:
                    if FileFormat == 2 or FileFormat == 4:
                        interrupt = 0
                        try:
                            for i, _qmax in enumerate(array_qmax):
                                print(f'qmax = {_qmax}, loop {i} / {len(array_qmax)-1}' )
                                if FileFormat == 2:
                                    _thmax = np.arcsin(_qmax * lambd / 4/np.pi)*180*2/np.pi
                                    th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,_thmax,QStep,incoherent,elements,ratio) #diff is just needed to define the values of S at low q values
                                    Imeas = np.array(Imeas)
                                    Ibkg = np.array(Ibkg)
                                    Ibrut=Imeas-bkfactor*Ibkg
                                
                                elif FileFormat == 4:
                                    q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,_qmax,incoherent,elements,ratio) 
                                    Imeas = np.array(Imeas)
                                    Ibkg = np.array(Ibkg)
                                    Ibrut=Imeas-bkfactor*Ibkg
                                
                                R=np.arange(0,rmax,0.0005)
                                R[0]=0.000001
                                q_old=q
                                pasq = q_old[1]-q_old[0]   
                                try:
                                    addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                                except: raise ValueError(('/\\'*30+"\n\nError, check q range\n\n"+'/\\'*30))
                                q=np.concatenate((addq,q_old))
                                
                                S = np.zeros((eq.prout+1,len(q)))
                                gr = np.zeros((eq.prout+1,len(R)))
                                F = np.zeros((eq.prout+1,len(R)))
                                delF = np.zeros((len(R)))
                                Sinv = np.zeros((len(q)))
                                
                                diff = np.abs(q-q_old[0]).argmin()
                                
                                ind = (np.abs(R - rmin)).argmin()
                                now = time.time()
                                args = (q,S,F,gr,Imeas,Ibkg,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                                out = minimize(eq.chisquare_DAC_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})

                                out_rhomi = out.params["rhomi"].value
                                out_bkfactor = out.params["bkfactor"].value
                                chi2function = out.residual[0]
                                print(f'qmax = {_qmax:.3f} rhomi = {out_rhomi:.3f} bkfactor = {out_bkfactor:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                                results[i] = [_qmax, out_rhomi, out_bkfactor, chi2function]
                                if debug: print(fit_report(out))
                                interrupt = interrupt +1
                        except KeyboardInterrupt:
                            results = results[0:interrupt,:]
                            print()
                            print('Interrupted !')
                            print()
                                    
                        
                    elif FileFormat == 3 or FileFormat == 1:
                        interrupt = 0
                        try:              
                            for i, _qmax in enumerate(array_qmax):
                                print(f'qmax = {_qmax}, loop {i} / {len(array_qmax)-1}' )
                                if FileFormat == 3:
                                    q,Ibrut, incohP, f2med, fmed2=eq.openFiles_3(filename,qmin,_qmax,incoherent,elements,ratio) 
                                    
                                elif FileFormat == 1:
                                    _thmax = np.arcsin(_qmax * lambd / 4/np.pi)*180*2/np.pi
                                    th,q,Ibrut, incohP, f2med, fmed2=eq.openFiles_1(filename,lambd,thmin,_thmax,QStep,incoherent,elements,ratio)
                        
                                R=np.arange(0,rmax,0.0005)
                                R[0]=0.000001
                                q_old=q
                                pasq = q_old[1]-q_old[0]
                                try:
                                    addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1]
                                except: raise ValueError(('/\\'*30+"\n\nError, check q range by visualize [v]\n\n"+'/\\'*30))
                                q=np.concatenate((addq,q_old))
                                
                                S = np.zeros((eq.prout+1,len(q)))
                                gr = np.zeros((eq.prout+1,len(R)))
                                F = np.zeros((eq.prout+1,len(R)))
                                delF = np.zeros((len(R)))
                                Sinv = np.zeros((len(q)))
                                
                                diff = np.abs(q-q_old[0]).argmin()
                                
                                ind = (np.abs(R - rmin)).argmin()
                                now = time.time()
                                args = (q,S,F,gr,Ibrut,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                                out = minimize(eq.chisquare_LV_LMFIT, fit_params, args = args, method='nelder', options = {'maxfev':500,'fatol':0.0001,'xatol':0.0001})
                                out_rhomi = out.params["rhomi"].value
                                chi2function = out.residual[0]
                                if FileFormat == 1 or FileFormat == 3:
                                    print(f'qmax = {_qmax} rhomi = {out_rhomi:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )                        
                                if FileFormat == 2 or FileFormat == 4:    
                                    print(f'qmax = {_qmax} rhomi = {out_rhomi:.3f} bkfactor = {out_bkfactor:.3f} chi2 = {chi2function:.6f} Nelder-Mead took {(time.time()-now)/60:.5f} min' )
                                results[i] = [_qmax, out_rhomi, 0,chi2function]
                                if debug: print(fit_report(out))
                                interrupt = interrupt +1
                        except KeyboardInterrupt:
                            results = results[0:interrupt,:]
                            print()
                            print('Interrupted !')
                            print()               
                        
                    allPar = [eval(i) for i in ParListRef]
                    header = eq.header_x(allPar , qmax_min, qmax_max, initialinitialrhomi,FileFormat_names[FileFormat])
                    

                    """ATTENTION: if change name in savetxt change name in eq.rename_F45 too !!!!"""
                    
                    np.savetxt(f'{F5_save_filename}_{filename[0:-4]}-{filename_bkg[0:-3]}_{qmax_min}-{qmax_max}.txt',results, header = header)
                    print(f'The file with the results has been saved as: {F5_save_filename}_{filename[0:-4]}-{filename_bkg}_{qmax_min}-{qmax_max}.txt')
            
                    
                    if FileFormat == 2 or FileFormat == 4:
                        fig_title =  'Loop over qmax'
                        fig, ax1 = plt.subplots()
                        ax1.set_title(fig_title)
                        color = 'tab:red'
                        ax1.set_xlabel('qmax (nm$^{-1}$)')
                        ax1.set_ylabel('bkfactor', color=color)
                        ax1.plot(np.transpose(results)[0], np.transpose(results)[2], color=color)
                        ax1.tick_params(axis='y', labelcolor=color)
                        ax2 = ax1.twinx()
                        color = 'tab:blue'
                        ax2.set_ylabel('chisquare', color=color)
                        ax2.plot(np.transpose(results)[0], np.transpose(results)[3], color=color)
                        ax2.tick_params(axis='y', labelcolor=color)
                        fig.savefig(F5_save_filename+'bkfactor'+str(qmax_min)+'_'+str(qmax_max)+'.png',dpi=DPI)        

                        plt.title(fig_title)
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                    
                    fig_title =  'Loop over qmax'
                    fig, ax1 = plt.subplots()
                    ax1.set_title(fig_title)
                    color = 'tab:red'
                    ax1.set_xlabel('qmax (nm$^{-1}$)')
                    ax1.set_ylabel('rho', color=color)
                    ax1.plot(np.transpose(results)[0], np.transpose(results)[1], color=color)
                    ax1.tick_params(axis='y', labelcolor=color)
                    ax2 = ax1.twinx()
                    color = 'tab:blue'
                    ax2.set_ylabel('chisquare', color=color)
                    ax2.plot(np.transpose(results)[0], np.transpose(results)[3], color=color)
                    ax2.tick_params(axis='y', labelcolor=color)
                    fig.savefig(F4_save_filename+'_chisquare_'+str(qmax_min)+'_'+str(qmax_max)+'.png',dpi=DPI)
                    plt.title(fig_title)
                    plt.show()
                    plt.pause(0.001) #used to update figure in terminal mode
                    
                except Exception as e: print(e)
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)
            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print()  
            
    
    ############################################################         
            
    
        elif case in ('Chi2','c'):
            try:
                print('\nChi2 selected')
                print()

                if FileFormat == 2 or FileFormat == 4:
                    rmax = 0.4
                    if FileFormat == 2:
                        thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                        thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi
                        th,q,Imeas,Ibkg,incohP, f2med, fmed2=eq.openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio) #diff is just needed to define the values of S at low q values
                        Imeas = np.array(Imeas)
                        Ibkg = np.array(Ibkg)
                        Ibrut=Imeas-bkfactor*Ibkg
                    
                    elif FileFormat == 4:
                        q,Imeas,Ibkg, incohP, f2med, fmed2=eq.openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio) 
                        Imeas = np.array(Imeas)
                        Ibkg = np.array(Ibkg)
                        Ibrut=Imeas-bkfactor*Ibkg
        
                    
                    R=np.arange(0,rmax,0.0005) 
                    R[0]=0.000001
                    q_old=q
                    pasq = q_old[1]-q_old[0]   
                    addq=np.sort(np.arange(q_old[0],0,-pasq))[:-1] 
                    q=np.concatenate((addq,q_old))
                    
                    S = np.zeros((eq.prout+1,len(q)))
                    gr = np.zeros((eq.prout+1,len(R)))
                    F = np.zeros((eq.prout+1,len(R)))
                    delF = np.zeros((len(R)))
                    Sinv = np.zeros((len(q)))
                    
                    print('Select rmin')
                
                    try: 
                        rmin = float(input())
                        print('\nrmin = '+ str(rmin))
                        ind = (np.abs(R - rmin)).argmin()
                        diff = np.abs(q-q_old[0]).argmin()
                        
                        print('\nSelect rhomi center')
                        initrhomi = float(input())
                        print('\nSelect rhomi range (+-)')
                        rhomi_range = float(input())
                        
                        print('\nSelect bkfactor center')
                        initbkfactor = float(input())
                        print('\nSelect bkfactor range (+-)')
                        bkfactor_range = float(input())
                        
                        print('\nSelect pixels number')
                        pixels = float(input())
                        print()
                        
                        
                        x = np.linspace(initrhomi-rhomi_range/2,initrhomi+rhomi_range/2,int(np.sqrt(pixels)))
                        y = np.linspace(initbkfactor-bkfactor_range/2,initbkfactor+bkfactor_range/2,int(np.sqrt(pixels)))
                        
                        Z = np.zeros((len(x),len(y)))
                        tot = len(x)*len(y)
                        
                        print()
                        print(f'Real size = {len(x)} x {len(y)} = {tot} pixels')
                        print()
                        print('Starting @ '+time.strftime("%Y:%m:%d %H:%M:%S" ))
                        print()
                        fit_params = Parameters()
                        fit_params.add('rhomi', value = initrhomi, min = max(1,initrhomi-50), max = initrhomi+50)
                        fit_params.add('bkfactor', value = initbkfactor, min = 0.7, max = 1.2)
                        
                        args = (q,S,F,gr,Imeas,Ibkg,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv)
                        
                        print()
                        print('Control+C to interrupt')
                        print()
                        
        
                        try:
                            for i,xval in enumerate(x):
                                for j,yval in enumerate(y):
                                    fit_params['rhomi'].set(xval)
                                    fit_params['bkfactor'].set(yval)
                                    Z[i][j]=eq.chisquare_DAC_LMFIT(fit_params,*args)
                                    k = i*len(y)+j
                                    jj = (k + 1)
                                    eq.update_progress(jj,tot)
                        except KeyboardInterrupt:
                            print()
                            print('Interrupted !')
                            print()
                    
                        print()
                        print('Ended @ '+time.strftime("%Y:%m:%d %H:%M:%S" ))
                        print()
                        
                        extent=[x[0],x[-1], y[0],y[-1]]
                        levels = list(np.linspace(Z.min()*1.5,Z.max(),5))
                        fig_title =  'Chisquare'
                        plt.figure(fig_title)
                        plt.title(fig_title)
                        plt.imshow(Z.T, cmap = 'gray',origin='lower', extent=extent , aspect="auto")
                        plt.colorbar()
                        cset = plt.contour(Z.T,levels, linewidths=2,extent=extent, aspect='auto',colors = 'yellow',origin='lower')
                        plt.clabel(cset, inline=10, fontsize=20)
                        plt.xlabel('density (at/nm$^3$)')
                        plt.ylabel('bkfactor')
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
                        
                        norm = Z.T/Z.min()
                        
                        fig_title =  'Normalized Chisquare'
                        plt.figure(fig_title)
                        plt.title(fig_title)
                        im = plt.imshow(norm, cmap='seismic',origin='lower',extent=extent, aspect='auto')
                        plt.colorbar(im)
                        cset = plt.contour(norm,[1,5,10,20,50,100], linewidths=2,extent=extent, aspect='auto',colors = 'yellow')
                        plt.clabel(cset, inline=10, fontsize=20)
                        plt.xlabel('rhomi (at/nm$^3$)')
                        plt.ylabel('bkfactor')
                        plt.plot([extent[0],extent[1]],[initbkfactor,initbkfactor],linewidth=2, color='red')
                        plt.plot([initrhomi,initrhomi],[extent[2],extent[3]],linewidth=2, color='red')
                        plt.xlim(extent[0],extent[1])
                        plt.ylim(extent[2],extent[3])
                        plt.show()
                        plt.pause(0.001) #used to update figure in terminal mode
            
                        print()
                        print('\nrmin = '+ str(rmin))
                        print()
                        print('\nrhomi = '+ str(initrhomi))
                        print()
                        print('\nbkfactor = '+ str(initbkfactor))
                        print()
                        np.savetxt(F3_save_filename+'_Chi2.txt',Z)
                        np.savetxt(F3_save_filename+'_Chi2_norm.txt',norm)
                    except Exception as e: print(e)
                else:
                    print()
                    print('Chi2 available only with background, FileFormat = 2 or 4')
                    print()
            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print()  
            except Exception as e: print(e)
            
    ############################################################    
                
        elif case in ('Par','Parameters','P','p','par'):
            try:
                while True:
                    print('/\\'*30)
                    print()
                    print(f'[01] {ParListNames[0]} = '+filename)
                    if FileFormat == 2 or FileFormat == 4:  
                        print(f'[02] {ParListNames[1]} = '+filename_bkg)
                    print(f'[03] {ParListNames[2]} = {FileFormat_names[int(FileFormat)]}')
                    if debug or FileFormat == 2 or FileFormat == 4:  
                        print(f'[04] {ParListNames[3]} = {incoherent}')
                    if FileFormat == 1 or FileFormat == 2:
                        print(f'[05] {ParListNames[4]}  = {lambd}')
                    print(f'[06] {ParListNames[5]} = {qmin}')
                    print(f'[07] {ParListNames[6]} = {qmax}')
                    print(f'[08] {ParListNames[7]} = '+ str(rmin))
                    if FileFormat == 2 or FileFormat == 4:  
                        print(f'[09] {ParListNames[8]}  = '+ str(bkfactor))
                    print(f'[10] {ParListNames[9]} = '+ str(rhomi))
                    print(f'[11] {ParListNames[10]} = {Normalization}')
                    print(f'[12] {ParListNames[11]} = {SofQ }')
                    print(f'[14] {ParListNames[13]} = {F3_save_filename}')
                    print(f'[15] {ParListNames[14]} = {F4_save_filename}')
                    print(f'[16] {ParListNames[15]} = {F5_save_filename}')
                    print(f'[17] {ParListNames[16]} = {elements}')
                    print(f'[18] {ParListNames[17]} = {ratio}')
                    print()
                    
                    print('Select parameter number to change, i to create Init.txt, q to exit')
                    case = input()
                    print()
                    if case in ('i','I','init','Init'):
                        print()
                        allPar = [eval(i) for i in ParListRef]
                        eq.save_init_file('Init.txt',ParListNames,allPar)
                        print()
                    elif case == 'q':
                        break #exit Parameters loop and go back to main loop
                    else:
                        try: 
                            case = int(case)-1
                            if  not 0 <= case <= 17: raise ValueError('goto except')
                            
                            print(f'parameter = {case+1}, name = {ParListNames[case]}, value = {eval(ParListRef[case])}')
                            if case == 2:
                                print(f"""
                1 => {FileFormat_names[1]}
                2 => {FileFormat_names[2]}
                3 => {FileFormat_names[3]}
                4 => {FileFormat_names[4]}
                """)
                            if case == 16:
                                print(elements)
                                elements,ratio = eq.change_el(elements,ratio)
                                print(elements)
                            elif case == 17:
                                elements,ratio= eq.change_ratio(elements,ratio)
                            elif case == 2: #210928 Section added, FileFormat must be integer
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = int(newval)
                                    if  not 1 <= newval <= 4: raise ValueError('goto except')
                                    exec(ParListRef[case]+" = " + str(newval))
                                    if debug == False:
                                        if FileFormat == 1 or FileFormat == 3:
                                            incoherent = 1
                                except:
                                    print('ERROR: FileFormat must be integer, try again !')
                                    print()
                            elif case == 3: #210928 Section added, incoherent must be integer = 0 or 1
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = int(newval)
                                    if  not 0 <= newval <= 1: raise ValueError('goto except')
                                    exec(ParListRef[case]+" = " + str(newval))
                                except:
                                    print('ERROR: incoherent must be 0 or 1, try again !')
                                    print()
                            elif case == 10: #210928 Section added, Normalization must be integer = 0 or 1
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = int(newval)
                                    if  not 0 <= newval <= 1: raise ValueError('goto except')
                                    exec(ParListRef[case]+" = " + str(newval))
                                except:
                                    print('ERROR: Normalization must be 0 or 1, try again !')
                                    print()
                            elif case == 11: #210928 Section added, SofQ must be integer = 0 or 1
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = int(newval)
                                    if  not 0 <= newval <= 1: raise ValueError('goto except')
                                    exec(ParListRef[case]+" = " + str(newval))
                                except:
                                    print('ERROR: SofQ must be 0 or 1, try again !')
                                    print()
                            elif case == 12: #210928 Section added, QStep must be integer = 0 or 1
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = int(newval)
                                    if  not 0 <= newval <= 1: raise ValueError('goto except')
                                    exec(ParListRef[case]+" = " + str(newval))
                                except:
                                    print('ERROR: QStep must be 0 or 1, try again !\nUsed only for debugging purpose')
                                    print()
                            elif 4 <= case <= 9:
                                print('Enter new numerical value')
                                newval = input()
                                print()
                                try:
                                    newval = float(newval)
                                    exec(ParListRef[case]+" = " + str(newval))
                                except:
                                    print('ERROR: Not a number, try again !')
                                    print()
                            else:
                                print('Enter new value')
                                newval = input()
                                print()
                                exec(ParListRef[case]+" = " + '"' + newval +'"')
                                if case == 13:
                                    F3_save_filename_ref = F3_save_filename
                                    F3_save_pointer = int(0)
                                if case == 14:
                                    F4_save_filename_ref = F4_save_filename
                                    F4_save_pointer = int(0)
                                if case == 15:
                                    F5_save_filename_ref = F5_save_filename
                                    F5_save_pointer = int(0)

                                    
                            print(f'parameter = {case+1}, name = {ParListNames[case]}, value = {eval(ParListRef[case])}')
                            print()
                        except:
                            print()
                            print('ERROR: This is not a valid selection, try again !')
                            print()
                
                    
                thmin= np.arcsin(qmin * lambd / 4/np.pi)*180*2/np.pi
                thmax= np.arcsin(qmax * lambd / 4/np.pi)*180*2/np.pi

            except KeyboardInterrupt:
                    print()
                    print('Interrupted !')
                    print()  
            except Exception as e:
                if debug:
                    print(traceback.format_exc())
                else:
                    print(e)

    ############################################################ 
    
        elif case in ('About','a','A'):
            print('/\\'*30)
            print("""
Amorpheus: a python-based software for the treatment of x-ray scattering data of amorphous and liquid systems
Version 1.0.5

Contacts: Silvia Boccato (silvia.boccato@cnrs.fr), Yiuri Garino (yiuri.garino@cnrs.fr)

Copyright (c) 2020-2022 Silvia Boccato, Yiuri Garino, Guillaume Morard and Chrystele Sanloup

Download: https://github.com/CelluleProjet/Amorpheus
How to cite: https://doi.org/10.1080/08957959.2022.2032032

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 724690).
""")
            print('/\\'*30)
            print()

    ############################################################ 
    
        elif case in ('Quit','q','Q'):
            print('/\\'*30)
            print()
            print('\nSee you again !!!')
            print()
            print('Quitting @ ' + time.strftime("%Y/%m/%d %H:%M:%S"))
            print()
            break
    
    ############################################################ 
    
        elif case in ('m'):
            local_vars = list(locals().items())
            print()
            print('Variable and size in bytes')
            print()
            mem = 0
            for var, obj in local_vars:
                print(var, sys.getsizeof(obj))
                mem = mem +  sys.getsizeof(obj)
            print()
            print(f'{len(local_vars)} Variables for {eq.convert_size(mem)}')
            print()
    
            
        else:
            print('\nNo matching Case, try again')
    
    except Exception as e:
        if debug:
            print(traceback.format_exc())
        else:
            print(e)
    except KeyboardInterrupt:
            print()
            print('The only way to exit the Continuum is Q')
            print()
