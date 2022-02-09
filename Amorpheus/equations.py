# -*- coding: utf-8 -*-
"""

Equations for Amorpheus
Version 1.0.5

Contacts: Silvia Boccato (silvia.boccato@cnrs.fr), Yiuri Garino (yiuri.garino@cnrs.fr)

Copyright (c) 2020-2022 Silvia Boccato, Yiuri Garino, Guillaume Morard and Chrystele Sanloup

Download: https://github.com/CelluleProjet/Amorpheus

This file is part of Amorpheus.
Amorpheus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Amorpheus is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Amorpheus.  If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 724690).

"""
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import time
import configparser
import ast
import math
from ezodf import opendoc, Sheet

from scipy.interpolate import Akima1DInterpolator
from numpy import linspace
import os

import xraylib as xr
import periodictable 

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

prout=4

symbol = ['♠','♣', '♥', '♦', '←', '↑','→', '↓', '↔', '«', '»', '≡', '¡', '¿', '¤','¶' ,'‡', '†', '§','◊','∫']

def update_progress(run,total, barLength = 20):
    progress = run/total
    block = int(round(barLength*progress))
    text = "\r"+time.strftime("%Y:%m:%d %H:%M:%S") + " [{0}] {1}%".format( symbol[block]*block + " "*(barLength-block), int(progress*100)) +  f' Run {run} over {total}'
    sys.stdout.write(text)
    sys.stdout.flush()
    
def getSandF(q,S,F,gr,Ibrut,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ, delF, Sinv):        
    S[0][diff:]=Normalize(Ibrut,q[diff:],rhomi,ratio,incohP,f2med, fmed2)
    
    if SofQ == 0:
        S[0][0:diff]=0
    elif SofQ == 1:
        S[0][0:diff]=S[0][diff]
    elif SofQ == 2:
        S[0][0:diff]=prolong(q,S[0],diff,1) 
    
    F[0] = calcF_speed(R,q,S[0])   
    gr[0]=1+F[0]/(4*np.pi*R*rhomi)
    
    ind = (np.abs(R - rmin)).argmin()
    delF = np.array(F[0]+4*np.pi*R*rhomi)
    for i in range(0,prout):   
        Sinv = 1/q*invdelF_speed(q,R,delF,ind) 
        S[i+1] = S[i]*(1-Sinv) 
        F[i+1] = calcF_speed(R,q,S[i+1])
        delF=F[i+1]+4*np.pi*R*rhomi
        gr[i+1]=1+F[i+1]/(4*np.pi*R*rhomi)
    return S[i+1],F[i+1],gr[i+1]

def lorenzianasenzaretta(x, a, b, gamma, x0):
    return a+b/(3.14*gamma*(1+((x-x0)/gamma)**2))

def prolong(q,S0,diff,plotprolong=None):
    p0=[80,20,.1,30]
    proprova=850
    ac, pcovUS = curve_fit(lorenzianasenzaretta,q[diff:diff+proprova],S0[diff:diff+proprova],p0,maxfev=30000)
    if plotprolong==1:
        plt.plot(q[0:diff+proprova],S0[0:diff+proprova])
        plt.plot(q,lorenzianasenzaretta(q,ac[0],ac[1],ac[2],ac[3]))
        plt.plot(q[0:diff],lorenzianasenzaretta(q,ac[0],ac[1],ac[2],ac[3])[0:diff])
        plt.show()
        print(f'q(0)={q[0]}, S0(0)={lorenzianasenzaretta(q,ac[0],ac[1],ac[2],ac[3])[0]}')
    return lorenzianasenzaretta(q,ac[0],ac[1],ac[2],ac[3])[0:diff]
    

def getSandF3(q,S,F,gr,Ibrut,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,delF, Sinv):        
    S[0][diff:]=Normalize(Ibrut,q[diff:],rhomi,ratio,incohP,f2med, fmed2)
    
    if SofQ == 0:
        S[0][0:diff]=0
    elif SofQ == 1:
        S[0][0:diff]=S[0][diff]
    elif SofQ == 2:
        S[0][0:diff]=prolong(q,S[0],diff,1) 
    
    F[0] = calcF_speed(R,q,S[0])   
    gr[0]=1+F[0]/(4*np.pi*R*rhomi)
    
    ind = (np.abs(R - rmin)).argmin()
    delF = np.array(F[0]+4*np.pi*R*rhomi)
    for i in range(0,prout):   
        Sinv = 1/q*invdelF_speed(q,R,delF,ind) 
        S[i+1] = S[i]*(1-Sinv) 
        F[i+1] = calcF_speed(R,q,S[i+1])
        delF=F[i+1]+4*np.pi*R*rhomi
        gr[i+1]=1+F[i+1]/(4*np.pi*R*rhomi)
    return S,F,gr

def Normalize(Isa,q,rhomi,ratio,incohP,f2med, fmed2):
    f1int = np.trapz( (q/10)**2*(f2med+incohP)/fmed2 ,q) 
    f2int = np.trapz( Isa*(q/10)**2/fmed2 ,q)
    alpha = (-2*np.pi**2*rhomi*10**(-3)+f1int)/f2int 
    S = (alpha*Isa-incohP-(f2med-fmed2))/fmed2  #equation 10 from Morard 2013, HPR
    return S

def ScatF(q,incoherent,elements,ratio):
    Znumber = [ eval('periodictable.' + i + '.number') for i in elements]
    FF_Rayl = np.zeros((len(Znumber),len(q)))
    SF_Compt  = np.zeros((len(Znumber),len(q)))
    for line, Z in enumerate(Znumber):
        FF_Rayl[line] = np.array([xr.FF_Rayl(Z,float(i/4/np.pi/10)) for i in q])
        SF_Compt[line]  = np.array([xr.SF_Compt(Z,float(i/4/np.pi/10)) for i in q])
    if incoherent == 0:
        incohP = 0
    else:
        incohP= (ratio[:,None]*SF_Compt).sum(axis = 0)
    f2med = (ratio[:,None]*FF_Rayl**2).sum(axis = 0)
    fmed2 = ((ratio[:,None]*FF_Rayl).sum(axis = 0))**2
    return incohP, f2med, fmed2



def calcF_speed(r,Q,S):
    """ numpy broadcasting """
    return 2/np.pi*np.trapz(Q*(S-1)*(np.sin(r[:,None]*Q)),Q,axis=1)

def invdelF_speed(q,R,delF,ind):
    return np.trapz(delF[0:ind]*np.sin(R[0:ind]*q[:,None]),R[0:ind],axis=1)

def openFile(filename):
    """ ignore initial space and empty lines, select only number """
    with open(filename) as f:
        s = f.read().replace(',','.')
        lines = (line for line in StringIO(s) if line.strip() != "" and line.strip()[0].isdigit()) 
        FH = np.loadtxt(lines)
    
    return FH


def ModificationFunction(r,Q,modfa,modfb):
    deltaR=np.pi/Q[-1]*(1-np.exp(-np.abs(r[:,None]-modfa)/modfb))    
    return np.sin(Q*deltaR)/(Q*deltaR)
    
    
def calcF_speed_ModificationFunction(r,Q,S,modfa,modfb):
    ModF = ModificationFunction(r,Q,modfa,modfb)
    return 2/np.pi*np.trapz(ModF*Q*(S-1)*(np.sin(r[:,None]*Q)),Q,axis=1)  


def arange(start,stop,step):
    if (stop-start)%step == 0:
        n = int((stop-start)/step+1)
        _range = np.linspace(start,stop,n)
    else:
        _range = np.arange(start, stop, step)
    return _range


def test():
    print('q ( nm$^{-1}$ )')


def plotFile(filename, filename_bkg, fileformat, lambd, q_min, q_max):
    """ ignore initial space and empty lines, select only number """
    FH = openFile(filename)
    
    fig_title =  'Data from file'
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_title(fig_title)
    ax.plot(FH[:,0],FH[:,1],label = f'{filename}')
    if fileformat == 2 or fileformat == 4:
        FH_bk = openFile(filename_bkg)
        ax.plot(FH_bk[:,0],FH_bk[:,1],label = f'{filename_bkg}')
    if fileformat == 3 or fileformat == 4:
        plt.xlabel(r'q ( $\AA^{-1}$ )')
        functions=(lambda x: x * 10, lambda x: x / 10)
    else:
        plt.xlabel('2 $\Theta$ ( $\degree$ )') 
        functions=(lambda x : 4*np.pi*np.sin(np.pi/180*x/2)/lambd, lambda x : np.arcsin(x * lambd / 4/np.pi)*180*2/np.pi)
    plt.ylabel('Data from file')
    plt.grid()
    secax = ax.secondary_xaxis('top', functions=functions)
    secax.set_xlabel('q ( nm$^{-1}$ )')
    plt.legend()
    plt.show()
    plt.pause(0.001) 
    
    if fileformat == 1 or fileformat == 2:
        _q_min = 4*np.pi*np.sin(np.pi/180*FH[0,0]/2)/lambd
        _q_max = 4*np.pi*np.sin(np.pi/180*FH[-1,0]/2)/lambd
        
        while q_min < _q_min or q_max > _q_max: 
            if q_min < _q_min:
                print(f'\nError, qmin must be >= of {_q_min} nm\u207B\u00B9')
                print(f"\nEnter New Value, enter to accept '{_q_min}'")
                val = input()
                if not val:
                    q_min = _q_min
                else:
                    try:
                        q_min = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
            if q_max > _q_max:
                print(f'\nError, qmax must be <= of {_q_max} nm\u207B\u00B9')
                print(f"\nEnter New Value, enter to accept '{_q_max}'")
                val = input()
                if not val:
                    q_max = _q_max
                else:
                    try:
                        q_max = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()                          
    
    if fileformat == 3 or fileformat == 4:
        _q_min = FH[0,0]*10
        _q_max = FH[-1,0]*10
        while q_min < FH[0,0]*10 or q_max > FH[-1,0]*10:
            if q_min < FH[0,0]*10:
                print(f'\nError, qmin must be >= of {FH[0,0]*10} nm\u207B\u00B9')
                print(f"\nEnter New Value in nm\u207B\u00B9, enter to accept '{FH[0,0]*10}'")
                val = input()
                if not val:
                    q_min = FH[0,0]*10
                else:
                    try:
                        q_min = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
            if q_max > FH[-1,0]*10:
                print(f'\nError, qmax must be <= of {FH[-1,0]*10} nm\u207B\u00B9')
                print(f"\nEnter New Value in nm\u207B\u00B9, enter to accept '{FH[-1,0]*10}'")
                val = input()
                if not val:
                    q_max = FH[-1,0]*10
                else:
                    try:
                        q_max = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
    return q_min, q_max, _q_min, _q_max

def checkFile(filename, filename_bkg, fileformat, lambd, q_min, q_max):
    """ ignore initial space and empty lines, select only number """
    FH = openFile(filename)
        
    if fileformat == 1 or fileformat == 2:
        _q_min = 4*np.pi*np.sin(np.pi/180*FH[0,0]/2)/lambd
        _q_max = 4*np.pi*np.sin(np.pi/180*FH[-1,0]/2)/lambd
        
        while q_min < _q_min or q_max > _q_max: 
            if q_min < _q_min:
                print(f'\nError, qmin must be >= of {_q_min} nm\u207B\u00B9')
                print(f"\nEnter New Value, enter to accept '{_q_min}'")
                val = input()
                if not val:
                    q_min = _q_min
                else:
                    try:
                        q_min = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
            if q_max > _q_max:
                print(f'\nError, qmax must be <= of {_q_max} nm\u207B\u00B9')
                print(f"\nEnter New Value, enter to accept '{_q_max}'")
                val = input()
                if not val:
                    q_max = _q_max
                else:
                    try:
                        q_max = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()                          
    
    if fileformat == 3 or fileformat == 4:
        _q_min = FH[0,0]*10
        _q_max = FH[-1,0]*10
        while q_min < FH[0,0]*10 or q_max > FH[-1,0]*10:
            if q_min < FH[0,0]*10:
                print(f'\nError, qmin must be >= of {FH[0,0]*10} nm\u207B\u00B9')
                print(f"\nEnter New Value in nm\u207B\u00B9, enter to accept '{FH[0,0]*10}'")
                val = input()
                if not val:
                    q_min = FH[0,0]*10
                else:
                    try:
                        q_min = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
            if q_max > FH[-1,0]*10:
                print(f'\nError, qmax must be <= of {FH[-1,0]*10} nm\u207B\u00B9')
                print(f"\nEnter New Value in nm\u207B\u00B9, enter to accept '{FH[-1,0]*10}'")
                val = input()
                if not val:
                    q_max = FH[-1,0]*10
                else:
                    try:
                        q_max = float(val)
                    except:
                        print('\nNot a number, try again !')
                        sys.exit()
    return q_min, q_max, _q_min, _q_max


def openFiles_1(filename,lambd,thmin,thmax,QStep,incoherent,elements,ratio):

    data = openFile(filename)
    th = data[:,0]
    Ibrut = data[:,1]

    spl = Akima1DInterpolator(th,Ibrut)

    
    array_th=np.linspace(thmin,thmax,2500)
    array_Ibrut=spl(array_th)

    if QStep == 0:
        array_q=4*np.pi*np.sin(np.pi/180*array_th/2)/lambd 
    elif QStep == 1:

        q1=4*np.pi*np.sin(np.pi/180*array_th[0]/2)/lambd
        qm=4*np.pi*np.sin(np.pi/180*array_th[len(array_th)-1]/2)/lambd
        pas=(qm-q1)/len(array_th)

        array_q=np.zeros(len(array_th))
        array_q[0]=q1
        for i in range(1,len(array_th)):
            array_q[i]=q1+pas*i

    incohP, f2med, fmed2 = ScatF(array_q,incoherent,elements,ratio)
    return array_th,array_q,array_Ibrut, incohP, f2med, fmed2

def openFiles_3(filename,qmin,qmax,incoherent,elements,ratio):

    data = openFile(filename)
    q = data[:,0]
    Ibrut = data[:,1]
    spl = Akima1DInterpolator(q*10,Ibrut)
    
    array_q = np.linspace(qmin,qmax,2500)
    
    array_Ibrut=spl(array_q)    
    incohP, f2med, fmed2 = ScatF(array_q,incoherent,elements,ratio)
    return array_q,array_Ibrut, incohP, f2med, fmed2

def openFiles_4(filename,filename_bkg,qmin,qmax,incoherent,elements,ratio):

    data = openFile(filename)
    q = data[:,0]
    Imeas = data[:,1]
    
    
    spl = Akima1DInterpolator(q*10,Imeas)
    

    array_q=np.linspace(qmin,qmax,2500)
    array_Imeas=spl(array_q)

    data = openFile(filename_bkg)
    q = data[:,0]
    Ibkg = data[:,1]
    spl = Akima1DInterpolator(q*10,Ibkg)
    
    array_Ibkg=spl(array_q)    
    incohP, f2med, fmed2 = ScatF(array_q,incoherent,elements,ratio)
    
    return array_q,array_Imeas,array_Ibkg,incohP, f2med, fmed2

def openFiles_2(filename,filename_bkg,lambd,thmin,thmax,QStep,incoherent,elements,ratio):

    data = openFile(filename)
    th = data[:,0]
    Imeas = data[:,1]

    
    spl = Akima1DInterpolator(th,Imeas)

    
    array_th=np.linspace(thmin,thmax,2500)
    array_Imeas=spl(array_th)

    if QStep == 0:
        array_q=4*np.pi*np.sin(np.pi/180*array_th/2)/lambd 
    elif QStep == 1:

        q1=4*np.pi*np.sin(np.pi/180*array_th[0]/2)/lambd
        qm=4*np.pi*np.sin(np.pi/180*array_th[len(array_th)-1]/2)/lambd
        pas=(qm-q1)/len(array_th)

        array_q=np.zeros(len(array_th))
        array_q[0]=q1
        for i in range(1,len(array_th)):
            array_q[i]=q1+pas*i

    data = openFile(filename_bkg)
    th = data[:,0]
    Imeasbkg = data[:,1]
    
    spl = Akima1DInterpolator(th,Imeasbkg)

    array_Ibkg=spl(array_th)    

    incohP, f2med, fmed2 = ScatF(array_q,incoherent,elements,ratio)
    return array_th,array_q,array_Imeas,array_Ibkg,incohP, f2med, fmed2


def input_rmin():
    _not_ok = True 
    while _not_ok:
        try:
            print('\nEnter rmin_min in nm (suggested for iron 0.1, max 0.4)')
            rmin_min = float(input())
            if rmin_min > 0.4:
                raise ValueError("Error, max 0.4 !")
            print('\nEnter rmin_max in nm (suggested for iron 0.2, max 0.4)')
            rmin_max = float(input())
            if rmin_max > 0.4:
                raise ValueError("Error, max 0.4 !")
            if rmin_max < rmin_min:
                raise ValueError("Error, rmin_max < rmin_min !")
            print('\nEnter the step (nm), suggested value 0.001 or 0.01')
            steprmin = float(input())
            array_rmin=arange(rmin_min,rmin_max,steprmin)
            
            if array_rmin.max() != rmin_max:
                print(f'\nreal rmin_max = {array_rmin.max():.6f}')
            if array_rmin.max() > 0.4:
                raise ValueError("Error, max 0.4 !")
            else:
                print(f'\n {len(array_rmin)} loops, continue [y] ?')
                _ = input()
                if _ == 'y' or not _ : _not_ok = False
        except Exception as e: print(e)
    return rmin_min, rmin_max, steprmin, array_rmin

def chisquare_DAC_LMFIT(pars,q,S,F,gr,Imeas,Ibkg,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv):
    vals = pars.valuesdict()
    rhomi = vals['rhomi']
    bkfactor = vals['bkfactor']
    Isa=Imeas-bkfactor*Ibkg
    S,F,gr=getSandF(q,S,F,gr,Isa,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ, delF, Sinv)
    return integrate.trapz((F[0:ind]+4*np.pi*R[0:ind]*rhomi )**2 ,R[0:ind])


def chisquare_LV_LMFIT(pars,q,S,F,gr,Ibrut,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ,ind, delF, Sinv):
    vals = pars.valuesdict()
    rhomi = vals['rhomi']
    S,F,gr=getSandF(q,S,F,gr,Ibrut,rhomi,R,rmin,rmax,diff,incohP,ratio,f2med, fmed2,SofQ, delF, Sinv)
    return integrate.trapz((F[0:ind]+4*np.pi*R[0:ind]*rhomi )**2 ,R[0:ind])

def rename_F3(filename_F3, filename_ref, pointer):
    while True:
        if os.path.isfile(filename_F3+'_S.txt'):
            print(f"Filename '{filename_F3}' arleady used, change ? y/n")
            rename = input()
        else:
            break
        if rename == 'n':
            break
        if rename == 'y' or not rename:      
            print(f"Enter New Filename NO EXTENSION, enter to accept '{filename_ref+'_'+str(pointer)}'")
            rename = input()
            
            if not rename:
                filename_F3 = filename_ref+'_'+str(pointer)
                pointer += 1
            else:
                filename_F3 = rename
                pointer = int(0)
                filename_ref = filename_F3
    return filename_F3, filename_ref, pointer

def rename_F45(filename_F45, filename_ref, pointer, filename, filename_bkg, ref1, ref2):
    while True:
        if os.path.isfile(f'{filename_F45}_{filename[0:-4]}-{filename_bkg[0:-3]}_{ref1}-{ref2}.txt'):
            print(f"Filename '{filename_F45}' arleady used, change ? y/n")
            rename = input()
        else:
            break
        if rename == 'n':
            break
        if rename == 'y' or not rename:      
            print(f"Enter New Filename NO EXTENSION, enter to accept '{filename_ref+'_'+str(pointer)}'")
            rename = input()
            
            if not rename:
                filename_F45 = filename_ref +'_'+str(pointer)
                pointer += 1
            else:
                filename_F45 = rename
                pointer = int(0)
                filename_ref = filename_F45
    return filename_F45, filename_ref, pointer

def ratio100(elements,ratio):
    while np.round(ratio.sum(),4) != 1:
        print('/\\'*30)
        print()
        print(f"Attention Total content not = 1 (= {ratio.sum():.03f})")
        print()
        print("Change [c] Autofill [A] or quit [q] ? A/q")
        case = input()
        if case == 'c' or case == 'C':
            [print(f'{el} content: {ratio[index]:.03f}') for index, el in enumerate(elements)]
            for index, el in enumerate(elements):
                print()
                print(f'{el} content: {ratio[index]:.03f}, enter new value')
                val = input()
                try:
                    if not val:
                        ratio[index] = 0
                    else:
                        ratio[index] = val
                except:
                    print('Not a number, try again !')
        elif case == 'A' or case == 'a':
            [print(f'[{index}] {el} content: {ratio[index]:.03f}') for index, el in enumerate(elements)]
            print()
            print("Select element to autofill")
            el100 = input()
            try:
                newval = int(el100)
                tot = 0
                for i,val in enumerate(ratio):
                    if i != newval:
                        tot = tot + val
                ratio[newval] = 1 - tot        
            except:
                print('Not a number, try again !')
        else: 
            print('/\\'*30)
            print()
            print('\nSee you again !!!')
            print()
            print('Quitting @ ' + time.strftime("%Y/%m/%d %H:%M:%S"))
            print()
            sys.exit()
        print()
        print("New Elements content")
        print()
        [print(f'{el} content: {ratio[index]:.03f}') for index, el in enumerate(elements)]
        print()
        print('/\\'*30)
    return(elements,ratio)

def change_ratio(elements,ratio):
    [print(f'{el} content: {ratio[index]:.03f}') for index, el in enumerate(elements)]
    for index, el in enumerate(elements):
        print()
        print(f'{el} content: {ratio[index]:.03f}, enter new value')
        val = input()
        try:
            if not val:
                ratio[index] = ratio[index]
            else:
                ratio[index] = val
        except:
            print('Not a number, try again !')
    elements,ratio = ratio100(elements,ratio)
    return(elements,ratio)
    
def change_el(elements,ratio):
    for index, el in enumerate(elements):
        print(f'Element {el}, enter new element')
        val = input()
        if not val:
            val = elements[index]
        try:
            print(f"New element '" + val +f"', known as {eval('periodictable.' + val + '.name')} , Z = {eval('periodictable.' + val + '.number')}")
            print()
            elements[index] = val
        except:
            print('Element unknown')
    
    return(elements,ratio)

def header_r(P, rmin_min, rmin_max, initialinitialrhomi,FF_name):
    if P[3]:
        inco = 'yes'
    else:
        inco = 'no'
    a = f'\
    rmin(nm) density(at/nm^3) bkfactor chisquare\n\
    Parameters analysis: filenames {P[0]}, rmin from {rmin_min:.3f} to {rmin_max:.3f}, initial density {initialinitialrhomi},\
    q from {P[5]:.3f} to {P[6]:.3f}, incoherent '
    b = f'\
    , SofQ {P[11]}\n\
    Parameters experiment: lambda {P[4]:.3f}, FileFormat {FF_name}, QStep'
    if P[12]:
        qstep = 'linear'
    else:
        qstep = 'real'
    normalization = f' Normalization {P[10]}'
    return a+inco+b+qstep + normalization

def header_x(P, qmin_min, qmin_max, initialinitialrhomi,FF_name):
    if P[3]:
        inco = 'yes'
    else:
        inco = 'no'
    a = f'\
    rmin(nm) density(at/nm^3) bkfactor chisquare\n\
    Parameters analysis: filenames {P[0]}, qmin from {qmin_min:.3f} to {qmin_max:.3f}, initial density {initialinitialrhomi},\
    q from {P[5]:.3f} to {P[5]:.3f}, incoherent '
    b = f'\
    , SofQ {P[11]}\n\
    Parameters experiment: lambda {P[4]:.3f}, FileFormat {FF_name}, QStep'
    if P[12]:
        qstep = 'linear'
    else:
        qstep = 'real'
    normalization = f' Normalization {P[10]}'
    return a+inco+b+qstep+normalization

def test2():
    exec('globals()["thmin"]=3')

def save_init_file(filename,ParListNames,allPar):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.add_section('Parameters')
    for i,name in enumerate(ParListNames):
        print(name + ' = ' + str(allPar[i]))
        config['Parameters'][name] = str(allPar[i])
    with open(filename, 'w') as configfile:
        config.write(configfile)

def printPar(ParListNames,allPar):
    for i,name in enumerate(ParListNames):
        print(name + ' = ' + str(allPar[i]))



def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return "%s %s" % (s, size_name[i])

def batch_check():
    doc = opendoc('Batch.ods')
    sheet = doc.sheets[0]
    data = np.array([[cell.value for cell in row] for row in sheet.rows()])
    lines_n = len(data[:,0])

    for i,value in enumerate(data[:,0][::-1]):
        if value == None:
            data = np.delete(data,lines_n-i-1, axis = 0)
        try:
            if value[0] == '#':
                print()
                print(f'Found # at row {lines_n-i-1}, row skipped')
                print()    
                data = np.delete(data,lines_n-i-1, axis = 0)
        except:
            pass
    pressures = data[1:,0]
    temperatures = data[1:,1]
    folders = data[1:,2]
    batch_loop_n = len(pressures)
    return data, temperatures, pressures, folders, batch_loop_n

def batch_save(folder):
    doc = opendoc('Batch.ods')
    doc.saveas(folder)