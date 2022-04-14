# -*- coding: utf-8 -*-
"""

Remove Peaks and Smooth (RPS) - Support tool for the software Amorpheus to pre-treat the data
Version 1.0.1

Principal author: Silvia Boccato (silvia.boccato@upmc.fr ), Yiuri Garino (yiuri.garino@cnrs.fr)
Copyright (c) 2020-2022 Silvia Boccato, Yiuri Garino, Guillaume Morard and Chrystele Sanloup

Free download: https://github.com/CelluleProjet/Amorpheus
How to cite: https://doi.org/10.1080/08957959.2022.2032032 

This file is part of Amorpheus.
Amorpheus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Amorpheus is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Amorpheus.  If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 724690).

"""

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import LorentzianModel, GaussianModel, VoigtModel


"""
da aggiungere per compilare:
import os
os.environ['MATPLOTLIBDATA'] = 'H:\venv\install\Lib\site-packages\matplotlib\mpl-data'
import tkinter.filedialog
"""
bg = "snow3"#"DeepSkyBlue4"#"DeepSkyBlue4"#


padx= 2
pady= 2
profile = "Loretz"
poly = False
fix_base = True
debug = False
resize = True #Resize False not working

if poly:
    from lmfit.models import PolynomialModel
    degree = 4 #range 2 - 4
    background = PolynomialModel(degree,prefix='line_')
else:
    from lmfit.models import LinearModel
    background = LinearModel(prefix='line_')

peak = LorentzianModel(prefix='lor_')

def change_fit(event):
    global peak
    print("_-"*30)
    print()
    print(Combo.current(), Combo.get())
    print()
    print("_-"*30)
    if Combo.current() == 0:
        peak = LorentzianModel(prefix='lor_')
    elif Combo.current() == 1:
        peak = GaussianModel(prefix='gau_')
    elif Combo.current() == 2: 
        peak = VoigtModel(prefix='voi_')
    elif Combo.current() == 3: 
        print("""
    No Fit procedure,
    Erase will just remove the data in 
    the range X min X max
    Press "Accept" to erase the data
    """)

def debugcheck_1(debug):
    if debug:
        print("_-"*30)
        print("canvas 1")
        print("_-"*30)
        print()
        print(f'data_fig x {len(data_fig.get_xdata())}')
        print(f'data_fig y {len(data_fig.get_ydata())}')
        print(f'substr x {len(substr.get_xdata())}')
        print(f'substr y {len(substr.get_ydata())}')
        print(f'smooth_Fig x {len(smooth_Fig.get_xdata())}')
        print(f'smooth_Fig y {len(smooth_Fig.get_ydata())}')
        print("_-"*30)
        print()
        print("_-"*30)
        
def debugcheck_2(debug):
    if debug:
        print("_-"*30)
        print("canvas 2")
        print("_-"*30)
        print(f'fit x {len(fit.get_xdata())}')
        print(f'fit y {len(fit.get_ydata())}')
        print(f'zoom x {len(zoom.get_xdata())}')
        print(f'zoom y {len(zoom.get_ydata())}')
        print(f'sub2 x {len(substr2.get_xdata())}')
        print(f'sub2 y {len(substr2.get_ydata())}')
        print(f'data_y_substr_tmp {len(data_y_substr_tmp)}')
        print("_-"*30)
        print()
        print("_-"*30)
    

def ExitApplication():
    MsgBox = tk.messagebox.askquestion ('Quitting ...','Are you sure you want to quit ?',icon = 'warning')
    if MsgBox == 'yes':
        root.quit()     # stops mainloop
        root.destroy()
    else:
        tk.messagebox.showinfo('Return','Going back')
    
def onclick(event):
    global Xpos
    if debug: print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    
    Xpos.set(f'{event.xdata:0.3f}')

def setMin():
    global Xmin, Xmax, data_x_zoom, data_y_zoom
    Xmin.set(Xpos.get())
    global line_min, canvas1
    line_min.set_xdata([float(Xpos.get()),float(Xpos.get())])
    canvas1.draw()
    global canvas2, zoom, fit
    b = np.ma.masked_outside(data_x, float(Xmin.get()), float(Xmax.get())).mask #less
    data_x_zoom = data_x[~b]
    data_y_zoom = data_y_substr[~b]
    zoom.set_xdata(data_x_zoom)
    zoom.set_ydata(data_y_zoom)
    if debug: print(f'Zoom Min = {data_x_zoom.min()}')
    global ref_Fig2
    a = [float(Xmin.get()),float(Xmax.get())]
    ref_Fig2.set_xlim(a)
    b = [data_y_zoom.min()-abs(data_y_zoom.min()*0.1),data_y_zoom.max()+abs(data_y_zoom.max()*0.1)]
    ref_Fig2.set_ylim(b)
    canvas2.draw()
    
def setMax():
    global Xmax, data_x_zoom, data_y_zoom
    Xmax.set(Xpos.get())
    global line_max, canvas1
    line_max.set_xdata([float(Xpos.get()),float(Xpos.get())])
    canvas1.draw()
    global canvas2, zoom
    b = np.ma.masked_outside(data_x, float(Xmin.get()), float(Xmax.get())).mask #less
    data_x_zoom = data_x[~b]
    data_y_zoom = data_y_substr[~b]
    zoom.set_xdata(data_x_zoom)
    zoom.set_ydata(data_y_zoom)
    if debug: print(f'Zoom Max = {data_x_zoom.max()}')
    global ref_Fig2
    a = [float(Xmin.get()),float(Xmax.get())]
    ref_Fig2.set_xlim(a)
    b = [data_y_zoom.min()-abs(data_y_zoom.min()*0.1),data_y_zoom.max()+abs(data_y_zoom.max()*0.1)]
    ref_Fig2.set_ylim(b)
    canvas2.draw()

def Erase():
    About_Erase_text = """
    No Fit procedure,
    Erase will just remove the data in 
    the range X min X max
    Press "Accept" to erase the data
    """
    tk.messagebox.showinfo("About_Erase", About_Erase_text)
        
def Fit():
    global data_x_zoom, data_y_zoom, fit, ref_Fig2, out, data_y_substr_tmp, sub
    if Combo.current() == 3:
        """ Erase condition """
        Erase()
    elif Combo.current() == 4:
        """ Replace condition """
        slope = (data_y_zoom[0]-data_y_zoom[-1])/(data_x_zoom[0]-data_x_zoom[-1])
        intercept = data_y_zoom[0]-slope*data_x_zoom[0]
        idx_min = np.abs(data_x - float(Xmin.get())).argmin()
        idx_max = np.abs(data_x - float(Xmax.get())).argmin()
        data_y_substr_tmp[idx_min:idx_max] = intercept + slope*data_x[idx_min:idx_max]
        fit.set_xdata(data_x[idx_min:idx_max])
        fit.set_ydata(data_y_substr_tmp[idx_min:idx_max])

    else:
        """ Fit Procedure """
        
        slope = (data_y_zoom[0]-data_y_zoom[-1])/(data_x_zoom[0]-data_x_zoom[-1])
        intercept = data_y_zoom[0]-slope*data_x_zoom[0]
        if poly:
            pars = background.make_params(c0=intercept, c1=slope, c2=0, c3=0,c4=0,c5=0,c6=0,c7=0)
        else:
            pars = background.make_params( slope = slope,intercept = intercept)
        center = data_x_zoom[data_y_zoom.argmax()]
        amplitude = data_y_zoom.max() - (intercept+slope*center)
        sigma =1
        beta =1
        pars += peak.make_params(amplitude  = amplitude, center = center, sigma = sigma, beta = beta)
        model = peak + background
        if fix_base and not poly:
            pars['line_slope'].vary = False
        if Combo.current() == 2: 
            pars['voi_gamma'].set(vary=True, expr='')
            pars['voi_gamma'].vary = True
            pars['voi_gamma'].expr = None
            pars['voi_gamma'].value = 1
            for par in pars.values():
                print(par)
        out = model.fit(data_y_zoom, pars, x=data_x_zoom)
        print('/\\'*32)
        print()
        print(out.fit_report())
        print()
        print('/\\'*32)
        fit_x = np.linspace(data_x_zoom.min(),data_x_zoom.max(),1000)
        fit.set_xdata(fit_x)
        fit.set_ydata(out.eval(x = fit_x))
        debugcheck_2(debug)
        canvas2.draw()

        if fix_base:
            if poly:
                init_slope = background.eval(slope = slope,intercept = intercept, x = data_x)
                if degree == 2:
                    init_slope = background.eval(c0=intercept, c1=slope, c2=0, x = data_x)
                elif degree == 3:
                    init_slope = background.eval(c0=intercept, c1=slope, c2=0, c3=0, x = data_x)
                elif degree == 4:
                    init_slope = background.eval(c0=intercept, c1=slope, c2=0, c3=0,c4=0, x = data_x)
            else:
                init_slope = background.eval(slope = slope,intercept = intercept, x = data_x)
            fit_res = out.eval(x=data_x)
            res = np.zeros(len(fit_res), dtype=float)
            b = np.ma.masked_outside(data_x, float(Xmin.get()), float(Xmax.get())).mask
            res[~b] = data_y_substr[~b] - fit_res[~b] + init_slope[~b]
            res[b] = data_y_substr[b]
            data_y_substr_tmp = res 
        else:
            lor_ = out.eval_components(x=data_x).get('lor_') #data in all range
            lor__ = np.zeros(len(lor_), dtype=float) 
            b = np.ma.masked_outside(data_x, float(Xmin.get()), float(Xmax.get())).mask
            lor__[~b] = lor_[~b]
            data_y_substr_tmp = data_y_substr-lor__ #substraction only in fit range

    substr2.set_ydata(data_y_substr_tmp)
    debugcheck_2(debug)
    canvas2.draw()
    
def Show_Res():
    tk.messagebox.showinfo("Fit Results", out.fit_report())

def Accept():
    global data_y_substr, data_y_substr_tmp, sub,box_pts_input, data_y_smooth, data_y_save
    global canvas2, zoom, data_y_zoom
    global ref_Fig2
    global data_x, data_y
    if Combo.current() == 3:
        """ Erase condition """
        idx_min = np.abs(data_x - float(Xmin.get())).argmin()
        idx_max = np.abs(data_x - float(Xmax.get())).argmin()
        print(f'{idx_min} value = {data_x[idx_min]}')
        print(f'{idx_max} value = {data_x[idx_max]}')
        data_x = np.concatenate((data_x[0:idx_min],data_x[idx_max:-1]))
        data_y = np.concatenate((data_y[0:idx_min],data_y[idx_max:-1]))

        data_y_save = np.copy(data_y)
        data_x_zoom = np.copy(data_x)
        data_y_zoom = np.copy(data_y)
        data_y_substr = np.copy(data_y)
        data_y_smooth = np.copy(data_y)
        
        data_fig.set_xdata(data_x)
        data_fig.set_ydata(data_y)
        substr.set_xdata(data_x)
        substr.set_ydata(data_y)
        smooth_Fig.set_xdata(data_x)
        smooth_Fig.set_ydata(data_y)
        print(f'x min {data_x.min()} max {data_x.max()}')
        print(f'y min {data_y.min()} max {data_y.max()}')
        a = [data_x.min(),data_x.max()]
        ref_Fig1.set_xlim(a)
        b = [data_y.min(),data_y.max()]
        ref_Fig1.set_ylim(b)
        global line_min
        line_min.set_xdata([round(data_x.min(),3),round(data_x.min(),3)])
        Xmin.set(round(data_x.min(),3))
        
        
        zoom.set_xdata(data_x)
        zoom.set_ydata(data_y)
        fit.set_xdata(data_x)
        fit.set_ydata(data_y)
        substr2.set_xdata(data_x)
        substr2.set_ydata(data_y)
        ref_Fig2.set_xlim(ref_Fig1.get_xlim())
        ref_Fig2.set_ylim(ref_Fig1.get_ylim())
        global line_max
        line_max.set_xdata([round(data_x.max(),3),round(data_x.max(),3)])
        Xmax.set(round(data_x.max(),3))
        debugcheck_2(debug)
        canvas2.draw()
        debugcheck_1(debug)
        canvas1.draw()
    else:

        #resetting smooth on accept
        box_pts_input.set('1')
        data_y_smooth = smooth(data_y_substr,1)
        smooth_Fig.set_ydata(data_y_smooth)
        canvas1.draw()
        
        data_y_substr = data_y_substr_tmp
        substr.set_ydata(data_y_substr)
        debugcheck_1(debug)
        canvas1.draw()
        
        b = np.ma.masked_outside(data_x, float(Xmin.get()), float(Xmax.get())).mask #less
        data_x_zoom = data_x[~b]
        data_y_zoom = data_y_substr[~b]
        zoom.set_xdata(data_x_zoom)
        zoom.set_ydata(data_y_zoom)
        
        a = [float(Xmin.get()),float(Xmax.get())]
        ref_Fig2.set_xlim(a)
        b = [data_y_zoom.min()-abs(data_y_zoom.min()*0.1),data_y_zoom.max()+abs(data_y_zoom.max()*0.1)]
        ref_Fig2.set_ylim(b)
        fit.set_xdata(data_x_zoom)
        fit.set_ydata(data_y_zoom)
        debugcheck_2(debug)
        canvas2.draw()
        data_y_save = data_y_substr

def Reset():
    global data_y_substr, data_y, data_y_zoom
    data_y_substr = data_y
    substr.set_ydata(data_y_substr)
    canvas1.draw()
    substr2.set_ydata(data_y_substr)
    fit.set_xdata(data_x_zoom)
    fit.set_ydata(data_y_zoom)
    canvas2.draw()

def smooth(y, box_pts):
    """ 
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way 
    https://numpy.org/devdocs/reference/generated/numpy.convolve.html
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same') # same size then y, boundary effects
    #y_smooth_2 = np.convolve(y, box, mode='valid') # only full overlap, for debug purpose
    y_smooth[0:box_pts//2-1] = y[0:box_pts//2-1] # replacing left side boundary effects with original data
    y_smooth[-(box_pts-1-box_pts//2)+1:] = y[-(box_pts-1-box_pts//2)+1:] # replacing right side boundary effects with original data
    return y_smooth #y_smooth_2 # for debug purpose
    
def Smooth():
    """Button event"""
    global data_y_smooth, data_y_substr, data_y_save
    try: 
        box_pts_val = int(box_pts_input.get())
        print(box_pts_input.get())
        print(box_pts_val)
        if box_pts_val <= 0:
            tk.messagebox.showinfo("Smooth value error", "Please eneter a positive integer")
        if box_pts_val > 0:
            data_y_smooth = smooth(data_y_substr,box_pts_val)
            smooth_Fig.set_ydata(data_y_smooth)
            debugcheck_1(debug)
            canvas1.draw()
            data_y_save = data_y_smooth
    except:
       tk.messagebox.showinfo("Smooth value error2", "Unknown error")

def open_file():
    global data_x, data_y, data_x_zoom, data_y_zoom, data_y_substr, data_y_substr_tmp, data_y_smooth
    global data_fig, substr, smooth
    global zoom, fit, ref_Fig1, ref_Fig2
    
    filename = tk.filedialog.askopenfilename()
    print(filename)
    data_x = np.loadtxt(filename)[:,0]
    data_y = np.loadtxt(filename)[:,1]
    
    data_x_zoom = np.copy(data_x)
    data_y_zoom = np.copy(data_y)
    data_y_substr = np.copy(data_y)
    data_y_substr_tmp = np.copy(data_y)
    data_y_smooth = np.copy(data_y)
    
    data_fig.set_xdata(data_x)
    data_fig.set_ydata(data_y)
    substr.set_xdata(data_x)
    substr.set_ydata(data_y)
    smooth_Fig.set_xdata(data_x)
    smooth_Fig.set_ydata(data_y)
    print(f'x min {data_x.min()} max {data_x.max()}')
    print(f'y min {data_y.min()} max {data_y.max()}')
    a = [data_x.min(),data_x.max()]
    ref_Fig1.set_xlim(a)
    b = [data_y.min(),data_y.max()]
    ref_Fig1.set_ylim(b)
    global line_min
    line_min.set_xdata([round(data_x.min(),3),round(data_x.min(),3)])
    Xmin.set(round(data_x.min(),3))
    
    
    zoom.set_xdata(data_x)
    zoom.set_ydata(data_y)
    fit.set_xdata(data_x)
    fit.set_ydata(data_y)
    substr2.set_xdata(data_x)
    substr2.set_ydata(data_y)
    ref_Fig2.set_xlim(ref_Fig1.get_xlim())
    ref_Fig2.set_ylim(ref_Fig1.get_ylim())
    global line_max
    line_max.set_xdata([round(data_x.max(),3),round(data_x.max(),3)])
    Xmax.set(round(data_x.max(),3))
    print(f'Fig x lim {ref_Fig1.get_xlim()}')
    print(f'Fig y lim {ref_Fig1.get_ylim()}')
    canvas2.draw()
    canvas1.draw()
    
def About():
    About_text = """

Remove Peaks and Smooth (RPS) - Support tool for the software Amorpheus to pre-treat the data
Version 1.0.1

Principal author: Silvia Boccato (silvia.boccato@upmc.fr ), Yiuri Garino (yiuri.garino@cnrs.fr)

Copyright (c) 2020-2022 Silvia Boccato, Yiuri Garino, Guillaume Morard and Chrystele Sanloup

Free download: https://github.com/CelluleProjet/Amorpheus

How to cite: https://doi.org/10.1080/08957959.2022.2032032 

This file is part of Amorpheus.
Amorpheus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Amorpheus is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Amorpheus.  If not, see <https://www.gnu.org/licenses/>.

Acknowledgements:
This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 724690).

"""
    tk.messagebox.showinfo("About", About_text)
    
def save_file():
    global data_y_save
    filename = tk.filedialog.asksaveasfile(mode='w', defaultextension=".txt", title="Save As (default extension .txt)", filetypes=(("Text file", "*.txt"),("Data file", "*.dat") ))
    if filename is None: #cancel button
        return
    print(filename)
    np.savetxt(filename, np.c_[data_x,data_y_save])
    filename.close()

root = tk.Tk()
root.option_add('*Dialog.msg.font', 'sans 10 bold')

root.wm_title("RPS: Remove Peaks & Smooth")

if resize:
    mainframe = ttk.Frame(root, borderwidth=5, relief="sunken", width=200, height=100)
    
else:
    mainframe = tk.Frame(root, borderwidth=5, relief="sunken", width=200, height=100)
    mainframe.config(bg=bg)
    
"""
Creating fake data for init
"""
if poly:
    if degree ==2:
        pars = background.make_params(c0=10, c1=.1, c2=0)
    elif degree == 3:
        pars = background.make_params(c0=10, c1=.1, c2=0, c3=0)
    elif degree == 4:
        pars = background.make_params(c0=10, c1=.1, c2=0, c3=0,c4=0)
else:  
    pars = background.make_params( slope = .1,intercept = 10)
pars += peak.make_params(amplitude  = 1500, center = 50, sigma = 5, beta = 40)
model = peak + background
data_x = np.linspace(0,100,1000)
data_y = model.eval(pars,x= data_x)

np.savetxt('FakeData.txt', np.c_[data_x,data_y])

data_x_zoom = np.copy(data_x)
data_y_zoom = np.copy(data_y)
data_y_substr = np.copy(data_y)
data_y_substr_tmp = np.copy(data_y)
data_y_smooth = np.copy(data_y)
data_y_save = np.copy(data_y)

Fig1 = plt.figure()

ref_Fig1 = Fig1.add_subplot(111)

data_fig, = ref_Fig1.plot(data_x,data_y, 'k',label = 'Raw Data')

substr, = ref_Fig1.plot(data_x_zoom,data_y_substr,color='orange', label = 'Sub')
smooth_Fig, = ref_Fig1.plot(data_x,data_y_smooth,color='green', label = 'Smooth')
line_min = ref_Fig1.axvline(x=data_x.min(), color = 'blue', linestyle = ':', label = 'Zoom Min')
line_max = ref_Fig1.axvline(x=data_x.max(), color = 'blue', linestyle = '--', label = 'Zoom Max')
ref_Fig1.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, loc = 'upper right')
ref_Fig1.grid()
ref_Fig1.set_xlabel('X axis')
ref_Fig1.set_ylabel('Intensity (a.u.)')
ref_Fig1.set_title('Data')

Xpos = tk.StringVar()
Xpos.set('0.000')
Xmin = tk.StringVar()
Xmax = tk.StringVar()
box_pts_input = tk.StringVar()
box_pts_input.set('1')


Xmin.set(data_x.min())
Xmax.set(data_x.max())

canvas1 = FigureCanvasTkAgg(Fig1, master = mainframe)
canvas1.get_tk_widget().grid(row=0, column=0, columnspan = 6, rowspan=25, padx = padx, pady = pady)
canvas1.draw()

toolbar_frame1=tk.Frame(mainframe)
toolbar_frame1.grid(row=25, column=0, columnspan = 6, sticky='NESW', padx = padx, pady = pady)

toolbar1 = NavigationToolbar2Tk(canvas1,toolbar_frame1)
toolbar1.grid(row=26,column=0, sticky='NESW')

mainframe.pack()

Fig2 = plt.figure()
ref_Fig2 = Fig2.add_subplot(111)
zoom, = ref_Fig2.plot(data_x_zoom, data_y_zoom, 'ok-', label = 'Zoom')
fit, = ref_Fig2.plot(data_x_zoom, data_y_zoom, 'r-', label = 'fit')
substr2, = ref_Fig2.plot(data_x,data_y_substr, color='orange',label = 'Sub')
ref_Fig2.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, loc = 'upper right')
ref_Fig2.grid()
ref_Fig2.set_xlabel('X axis')
ref_Fig2.set_ylabel('Intensity (a.u.)')
ref_Fig2.set_title('Zoom Area')

canvas2 = FigureCanvasTkAgg(Fig2, master = mainframe)
canvas2.get_tk_widget().grid(row=0, column=6, columnspan = 6,rowspan=25, padx = padx, pady = pady)
canvas2.draw()

toolbar_frame2=tk.Frame(mainframe)
toolbar_frame2.grid(row=25, column=6, columnspan = 6, sticky='NESW', padx = padx, pady = pady)

toolbar2 = NavigationToolbar2Tk(canvas2,toolbar_frame2)
toolbar2.grid(row=26,column=6, sticky='NESW')

canvas1.mpl_connect('button_press_event', onclick)

"""Buttons sections """

label_Cursor = ttk.Label(mainframe, text="Cursor", foreground = 'black', font='sans 10 bold', anchor='center')
label_Cursor.grid(column=0, row=27, columnspan=1, padx = padx, pady = pady, sticky = 'NESW')

label_Xpos = ttk.Label(mainframe, textvariable=Xpos, font='sans 10 bold', anchor='center')
label_Xpos.grid(column=1, row=27, columnspan=1, padx = padx, pady = pady, sticky = 'NESW')# sticky='NESW'

button_set_xmin = tk.Button(master=mainframe, text="MIN", command=setMin,  bg = 'blue', fg = "white", font='sans 10 bold')
button_set_xmin.grid(column=2, row=27, columnspan=1, sticky = 'NESW', padx = padx, pady = pady)

label_Xmin = ttk.Label(mainframe, textvariable=Xmin, foreground = 'blue', font='sans 10 bold', anchor='center')
label_Xmin.grid(column=3, row=27, columnspan=1, padx = padx, pady = pady, sticky = 'NESW')

vlist = ["Lorentz", "Gauss", "Voigt", "Erase", "Replace"]
Combo = ttk.Combobox(mainframe, values = vlist, font='sans 10 bold')
Combo.set("Lorentz")
Combo.bind("<<ComboboxSelected>>", change_fit )
Combo.grid(column=4, row=27, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_smooth = tk.Button(master=mainframe, text="Smooth →", command=Smooth, font='sans 10 bold')
button_smooth.grid(column=0, row=28, columnspan=1, sticky = 'NESW', padx = padx, pady = pady)

button_smooth_value = tk.Entry(master=mainframe, textvariable=box_pts_input, justify='center', font='sans 10 bold')
button_smooth_value.grid(column=1, row=28, columnspan=1, sticky = 'NESW', padx = padx, pady = pady)

button_set_xmax = tk.Button(master=mainframe, text="MAX", command=setMax,  bg = 'blue', fg = "white", font='sans 10 bold')
button_set_xmax.grid(column=2, row=28, columnspan=1, sticky = 'NESW', padx = padx, pady = pady)

label_Xmax = ttk.Label(mainframe, textvariable=Xmax, foreground = 'blue', font='sans 10 bold', anchor='center')
label_Xmax.grid(column=3, row=28, columnspan=1, padx = padx, pady = pady, sticky = 'NESW')

#https://unicode-table.com/en/1F60E/

button_Fit = tk.Button(master=mainframe, text="FIT", command=Fit, font='sans 10 bold', bg = "gold")
button_Fit.grid(column=4, row=28, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_Show_Res = tk.Button(master=mainframe, text="Show Results", command=Show_Res, font='sans 10 bold')
button_Show_Res.grid(column=6, row=27, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_Substract = tk.Button(master=mainframe, text="Accept ✔ ", command=Accept, font='sans 10 bold')
button_Substract.grid(column=8, row=27, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_About = tk.Button(master=mainframe, text="About", command=About, fg = 'black', font='sans 10 bold')
button_About.grid(column=10, row=27, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_open = tk.Button(master=mainframe, text="Open", command=open_file, fg = 'white', font='sans 10 bold',bg="RoyalBlue2")
button_open.grid(row = 28, column=6, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_save = tk.Button(master=mainframe, text="Save", command=save_file, fg = 'white', font='sans 10 bold', bg = "forest green")
button_save.grid(row = 28, column=8, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)

button_quit = tk.Button(master=mainframe, text="Quit ✘", command=ExitApplication, bg = 'red3', fg = 'white', font='sans 10 bold')
button_quit.grid(row = 28, column=10, columnspan=2, sticky = 'NESW', padx = padx, pady = pady)


if resize:
    for x in range(12):
        tk.Grid.columnconfigure(mainframe, x, weight=1)
    
    for y in range(28):
        tk.Grid.rowconfigure(mainframe, y, weight=1)

root.mainloop()
