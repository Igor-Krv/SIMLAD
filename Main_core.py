import pandas
import tkinter as tk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from tkinter import ttk
from datetime import (datetime,date)
import json
from os.path import exists
import os.path
import FastEuler
import ctypes
try: # >= win 8.1
    ctypes.windll.shcore.SetProcessDpiAwareness(2)
except: # win 8.0 or less
    ctypes.windll.user32.SetProcessDPIAware()
import os
from pathlib import Path

from Widget_core import Widget
from Laser_core import Laser
from Driver_core import Driver
from Filter_core import Filter
from custom_pulse import Arbpulse
from json_core import config
from json_core import conf_pulsetrain
from tkinter import messagebox as tk_messagebox
from tkinter import filedialog
from ttkthemes import themed_tk
from ttkthemes import ThemedStyle

hbar = 1.054 *1e-25


debugtime = []
class App(ttk.Frame):
    def __init__(self, parent):
        ttk.Frame.__init__(self, parent)
        
        self.parent = parent 
        self.ready = 1          #non used #if there is no config file, all param are 0, can't calculate
        self.driver_px =1120
        self.driver_py = 60
        
        self.wdmx  = 530        #for wdm params on the main frame
        self.wdmy = 188
        self.pltsy = 10
        self.pltsx = 10
        self.slavestate = 0     #indicator if slave is activated
        
        self.datalstart = 9     #start and end of a laser data array 
        self.datalend = 38      # from 0 to 8 driver data is stored        
        self.intstart = 34      # start of interferometer data in Laser structure
        
        self.var = []           # array for calculation and plot options (both  master and slave )
        self.lbl = []
        self.myvar = []
        self.filter_driver =""  #filter type name
        self.dvar = []
            
        self.intvar =[]         #variables of interferometer for tkinter interface
        self.entry = []
        self.dentry = []        # entry interface for master wdm
        self.slavedentry = []   # entry interface for slave wdm
        self.intentry =[]
        self.slaveentry = []    # data for slave Laser Entry
        self.cdata =[]
           
        self.pltdata = []
        
        self.wdmvar = []        #variables of wdm for tkinter interface
        self.wdmdata = []       #1st wdm filter data
        self.setuplbl =[]       #lbl for wdm filter master
        
        self.slavewdmwar = []   #variables of slave wdm for tkinter interface
        self.slavewdmdata = []  #slave wdm filter data
        self.slavewdmlbl = []
        
        self.intdata = []
        self.intlbl =[]
        self.Imlbl = []         #array for images (setup scheme)

        self.photentry = []
        self.photlbl = []
        self.photvar = []
        self.photdata = np.zeros(3)
        
        self.slavedata = np.zeros(self.datalend)    #laser parameters of a slave
        self.data =  np.zeros(self.datalend)        #master laser parameters
        
        self.slavelbl = []
        self.slavevar = []
        
        
        
        self.mx = 0                         # variables for tracking mouse position 
        self.my = 0
        self.mx_p = 0
        self.my_p = 0
       
        self.datadstart = 0                 #driver array borders
        self.datadend = 9                   #from 9 to 33 laser data
        self.plt_settings = np.zeros(30)    #filetext for plot marks
        
        self.arbvar = []                    #arbitary pulses
        self.arbdata = np.zeros(self.datalstart) 
        self.arblbl =[] 
        self.arbentry =[]
        
        self.maxarbpulse = 100              #maximum number of arbitrary pulses
        self.arbtime = []
        self.arbfiltered = []
        self.pulseid = 0
        self.Arbpulses =[]
        self.arb_cache = np.zeros((10,self.maxarbpulse))
        
        self.slavepulseid = 0
        self.slaveArbpulses = []
        self.slavearb_cache =  np.zeros((10,self.maxarbpulse))
        
        
        self.transdata = np.zeros(2)
        
        self.responselabel = ttk.Label()
        self.response_wp_1 = ttk.Spinbox()   
        self.responsevar_1 = tk.StringVar()
        self.response_wp_2 = ttk.Spinbox()   
        self.responsevar_2 = tk.StringVar()

        self.stabdata = np.zeros(6)         #array for diagram of stability parameters
        for i in range(0,self.datalend):    #fill data and interface elemets arrays
            self.cdata.append(0)
            self.pltdata.append(0)
            self.wdmdata.append(0)
            self.slavewdmdata.append(0)
            self.intdata.append(0)
            
            self.myvar.append(tk.StringVar())       #variables for widgets
            self.wdmvar.append(tk.StringVar())
            self.slavewdmwar.append(tk.StringVar())
            self.intvar.append(tk.StringVar())
            self.photvar.append(tk.StringVar())
            
            self.entry.append(ttk.Spinbox())        #master laser array
            self.slaveentry.append(ttk.Spinbox())   #slave laser array
            self.dentry.append(ttk.Spinbox())       #master wdm array
            self.slavedentry.append(ttk.Spinbox())  #slave wdm array
            self.intentry.append(ttk.Spinbox())     #interferometer array
            self.photentry.append(ttk.Spinbox())    #photodetector array

             
            self.lbl.append(ttk.Label())            #labels for widgets
            self.setuplbl.append(ttk.Label())
            self.slavewdmlbl.append(ttk.Label())
            self.intlbl.append(ttk.Label())
            self.slavelbl.append(ttk.Label())
            self.Imlbl.append(ttk.Label())
            self.photlbl.append(ttk.Label())
                                
            self.slavevar.append(tk.StringVar())
            self.var.append(tk.BooleanVar())
            self.dvar.append(tk.BooleanVar())
            
            self.arbvar.append(tk.StringVar())
            self.arblbl.append(ttk.Label())
            self.arbentry.append(ttk.Spinbox())
        
        self.initUI()

        for i in range(0,self.maxarbpulse): #creating structures for arb pulses
            self.Arbpulses.append(Arbpulse(self.widgets.pulseframe_mod_aux))
            self.slaveArbpulses.append(Arbpulse(self.widgets.pulseframe_mod_aux)) 
         
            
        self.figure = plt.Figure(figsize=(16,9), dpi=100)       #plot display for main calculation
        self.canvas = FigureCanvasTkAgg(self.figure, self.widgets.pulseframe_reg)
        
        self.stabfigure = plt.Figure(figsize=(16,9), dpi=100)    #plot display for stability diagram
        self.stabcanvas = FigureCanvasTkAgg(self.stabfigure, self.widgets.pulseframe_stab)  
        
        self.modfigure = plt.Figure(figsize=(16,9), dpi=100)    #plot display for arbitrary pulses
        self.modcanvas = FigureCanvasTkAgg(self.modfigure, self.widgets.pulseframe_mod) 
    
        for i in range(0,self.datadend): #bind tk objects so values won't change unless you press <enter> to save them
                self.widgets.master_entry[i].spinboxm.bind("<FocusOut>", self.driver_reset)
                self.widgets.slave_entry[i].spinboxm.bind("<FocusOut>", self.driver_reset)
        for i in range(0, self.datalend):
            self.entry[i].bind("<FocusOut>", self.laser_reset)
            
            self.slaveentry[i].bind("<FocusOut>", self.laser_reset)
            
        for i in range(0, self.datalend):
            
            #self.slavedentry[i].bind("<FocusOut>", self.input_reset)
            self.intentry[i].bind("<FocusOut>", self.input_reset)
            self.photentry[i].bind("<FocusOut>", self.input_reset)
            self.slavedentry[i].bind("<FocusOut>", self.input_reset)
            self.dentry[i].bind("<FocusOut>", self.input_reset)
        for i in range(0,len(self.widgets.stabinput)):
            self.widgets.stabinput[i].spinboxm.bind("<FocusOut>", self.input_reset)
        for i in range(0,self.maxarbpulse):
            for j in range(0,10):
                self.Arbpulses[i].entry[j].bind("<FocusOut>", self.arbpulse_reset)
                self.slaveArbpulses[i].entry[j].bind("<FocusOut>", self.slavearbpulse_reset)
                
        self.option_changed()       #Initialise starting scheme image
        
    
    def maincalculation(self): #perform calculation in the "Regular Simulations" frame

        cm = [self.widgets.m_checkbox_clc[i].var.get() for i in range(0,self.widgets.master_clcoptions)]
        dm = [self.widgets.m_checkbox_plt[i].var.get() for i in range(0,self.widgets.master_pltoptions)]
        
        #cm - calculation master option, dm - display slave option
        
        # c0 = self.widgets.m_checkbox_clc[0].var.get()    #driver
        # c1 = self.widgets.m_checkbox_clc[1].var.get()    #laser output  N Q
        # c2 = self.widgets.m_checkbox_clc[2].var.get()    #phase
        # c3 = self.widgets.m_checkbox_clc[3].var.get()   #frequency 
        # c4 = self.widgets.m_checkbox_clc[4].var.get()    #temperature
        # c5 = self.widgets.m_checkbox_clc[5].var.get()    #noise
        # c6 = self.widgets.m_checkbox_clc[6].var.get()    #biexp temp
        # c7 = self.widgets.m_checkbox_clc[7].var.get()    #interference
        
        # d0 = self.widgets.m_checkbox_plt[0].var.get()   #driver
        # d1 = self.widgets.m_checkbox_plt[1].var.get()   #N Q
        # d2 = self.widgets.m_checkbox_plt[2].var.get()   #phase
        # d3 = self.widgets.m_checkbox_plt[3].var.get()   #frequency
        # d4 = self.widgets.m_checkbox_plt[4].var.get()   #temperature
        # d5 = self.widgets.m_checkbox_plt[5].var.get()   #biexp
        # d6 = self.widgets.m_checkbox_plt[6].var.get()   #interference
            
        cs = [self.widgets.s_checkbox_clc[i].var.get() for i in range(0,self.widgets.slave_clcoptions)]
        ds  = [self.widgets.s_checkbox_plt[i].var.get() for i in range(0,self.widgets.slave_pltoptions)]      
        #cs - calculation slave options, ds - display slave options
    
        #################################### master driver
        mIb = float(self.data[0])*1e-12
        mIp = float(self.data[1])*1e-12
        mrate = float(self.data[3])
        mbandwidth = float(self.data[4])
        mforder = int(self.data[5])
        mfp = float(self.data[6])
        mw = float(self.data[7])        
        mclctime = float(self.data[8])
        ##################################### slave driver
        slIb = float(self.slavedata[0])*1e-12
        slIp = float(self.slavedata[1])*1e-12
        slrate = mrate                        #slave has the same repetition rate
        slbandwidth = float(self.slavedata[4])
        slorder = int(self.slavedata[5])
        slfp = float(self.slavedata[6])
        slw = float(self.slavedata[7])        
        slctime = float(self.slavedata[8])
        slnpulses = float(self.slavedata[2])
        arg = np.zeros(len(self.data),dtype = np.float64)                            #parameters to pass in the master laser
        slavearg = []                         #parameters to pass in the slave laser
        for i in range(self.datalstart,self.datalend):
            arg[i - self.datalstart] = float(self.data[i]) 
        for i in range(self.datalstart,self.datalend):
            slavearg.append(float(self.slavedata[i]))
            
        arg[25] = float(self.intdata[0])   #interferometer delay
        arg[26] = float(self.intdata[1])   #theta 1
        arg[27] = float(self.intdata[2])   #theta 2
        slavearg[25] = arg[25]
        slavearg[26] = arg[26]
        slavearg[27] = arg[27]
        if ( self.data[2] != ""):           #checking pulse number
            npulses = float(self.data[2])

        if (npulses > 0):                                                    #calculation starts if there are at least 1 pulse                                      
            start = datetime.now()                                           #start time count 
            driver =  Driver(mIp,mIb,mfp,mw,mclctime,npulses,mrate)          #master driver
            print("master driver is ok")
            sldriver =  Driver(slIp,slIb,slfp,slw,slctime,slnpulses,slrate)  #slave driver
            wdm = Filter(mfp,mw,mclctime,npulses,mrate)                      #wdm both for the master and slave laser
            #signal from master driver
            sigtype = self.widgets.sigtype_menu.option_var.get()
            filtered = driver.sigfilter(driver,driver.pulses(driver,sigtype),mbandwidth,mforder,str(self.widgets.f_var.get()),sigtype )       
            
            if( slnpulses > 0):                                     
                slfiltered = sldriver.sigfilter(sldriver,sldriver.pulses(sldriver,sigtype),slbandwidth,slorder,str(self.widgets.f_var.get()),sigtype )
            #comparing time between slave and master (in case of difference in width and number)
            maxUppertime = driver.upperTime             #max time in seconds
            maxtime = driver.time                       #max time in points
            if(driver.upperTime < sldriver.upperTime):  #getting maximum time array from driver or slave
                maxtime = sldriver.time
                maxUppertime = sldriver.upperTime
            if(driver.upperTime > sldriver.upperTime):
                maxtime = driver.time
                maxUppertime = driver.time
            #print("Highest upper time " + str(maxUppertime) +" " + str(len(maxtime)) )            
            if (npulses < slnpulses ):                  #adding zeros to master's or slave's time depending on pulse number
                filtered = np.concatenate( [filtered, np.zeros(abs(len(filtered) - len(slfiltered)))*1e-2 ])
                driver.time = sldriver.time#np.concatenate( [driver.time, np.zeros(abs(len(driver.time) - len(sldriver.time)))*1e-2 ],axis = 0)
                driver.L = len(driver.time)
            elif (npulses > slnpulses):
                slfiltered = np.concatenate([slfiltered , np.zeros(abs(len(filtered) - len(slfiltered)))*1e-2 ])
                sldriver.time = driver.time#np.concatenate( [sldriver.time, np.zeros(abs(len(driver.time) - len(sldriver.time)))*1e-2 ],axis = 0)
                sldriver.L = len(sldriver.time)
            #TODO ! 
            laser = Laser(mIp,driver,arg)               #initiating lasers
            laserS = Laser(slIp,sldriver,slavearg)
            #print(f"slave laser time {len(laserS.time)}")
            debugtime = laserS.time
            output = []
            outputarr = []          #plot data
            outputlabels = []       #plot labels
            outputleg = []          #plot legends
            text_output = []        #temporary array to store output for saving
            name = ""

            self.widgets.threscurrM.var.set(round(laser.Ith*1e12,4)) #Threshold values
            self.widgets.threscurrS.var.set(round(laserS.Ith*1e12,4))
            if ( self.widgets.calc_var.get() == "Runge_Kutta"):
                
                if ( self.widgets.scipym.var.get()): #unstable, sci_ODE turned off
                    output = laser.sci_ODE(filtered,cm[2],cm[3],cm[5],cm[4],cm[7])
                    tempphase = output[2]
                    tempT = output[3]
                    nu = np.zeros(len(output[4]))
                    for i in range(0,len(output[4])-1):
                        nu[i] = (output[3][i+1] - output[3][i]) /(laser.Dt*2*np.pi)                                  
                else:
                    output = laser.runge_kutta(filtered,cm[2],cm[3],cm[5],cm[4],cm[7],cm[6])
                   
                print ("Runge Kutta scheme is chosen for slave master non-stochastic calculations")
            else:
                if (self.widgets.Euler_mod.option_var.get() == 'Regular'):
                    output = laser.euler_m(filtered,cm[2] , cm[3], cm[5] , cm[4],cm[7],0,0)
                    #outputE = FastEuler.fast_euler( driver.Ib, driver.Ip,filtered,driver.Dt,driver.L,driver.nPointsPerPulse,arg,cm[2],cm[3],cm[5],cm[4],cm[7],0,0)
                    #outputE[1] = np.sqrt(outputE[1]) * np.exp(1j*outputE[3])
                else: #unstable, turned off
                    output = FastEuler.fast_euler( driver.Ib, driver.Ip,filtered,driver.Dt,driver.L,driver.nPointsPerPulse,arg,cm[2],cm[3],cm[5],cm[4],cm[7],0,0)
                    output[1] = np.sqrt(output[1]) * np.exp(1j*output[3])
               
                print(f"Euler {self.widgets.Euler_mod.option_var.get()} scheme is chosen for slave master calculations")

            if(self.widgets.option_var.get().find('Master-WDM') != -1 ): 
                output[1] = wdm.wdm_sigfilter( wdm, output[1], float(self.wdmdata[0]), int(self.wdmdata[1]), float(self.wdmdata[2]),str(self.widgets.wdmvar1.get()))
            #photodecetor is measuring masters output
            # detector_outM = wdm.photodetector( wdm, output[1] , float(self.photdata[0]), int(self.photdata[1]), float(self.photdata[2]),"Bessel")                    
            #laser.E = output[1]#np.sqrt(output[1]) * np.exp(1j*output[3]) # in case method works with Q
            #np.sqrt(output[1])* np.exp(1j* output[3])
            detector_outM = wdm.photodetector( wdm, output[1]  , float(self.photdata[0]), int(self.photdata[1]), float(self.photdata[2]),"Bessel")             

            if (self.widgets.power_output.var.get()): #Power output in Watts
                Q = detector_outM
                #hbar = 1.054 *1e-25 #J*s
                #t = 0.001 #*1e-9 #s                 
                #lambd = 1550#*1e-9 #m #1550 #nm
                c = 2*1e17# m/s
                omega = c / laser.lambd*1e3
                #P = laser.eta  * laser.hbar * omega * Q / (2 * laser.Gamma * laser.tau_ph)
                detector_outM = ( laser.eta*laser.hbar *(c/(laser.lambd * 1e3) ) )/(2*laser.Gamma*laser.tau_ph)* Q                               
                
            if(self.slavestate == 1):
                if ( self.widgets.calc_var.get() == "Runge_Kutta"):
                    if ( self.widgets.scipym.var.get()): #unstable, turned off in Main core
                        outputs = laserS.sci_ODE_slave(slfiltered,output[1],cs[2],cs[3],cs[5],cs[4],cs[7],cs[6])
                        nuS = np.zeros(len(outputs[4]))
                        for i in range(0,len(outputs[3])-1):
                            nuS[i] = (outputs[3][i+1] - outputs[3][i]) /(laser.Dt*2*np.pi)               
                    else:
                        outputs = laserS.runge_kutta_slave(output[1], slfiltered,self.widgets.injfix.var.get(),cs[2],cs[3],cm[5],cm[4],cm[7],cm[6])
                else:
                    
                    if (self.widgets.Euler_mod.option_var.get() == 'Regular'):
                        outputs = laserS.euler_m_slave(output[1],sldriver,slfiltered,self.widgets.injfix.var.get(),cs[5],cs[4],cs[7],0,0,cs[2],cs[3],cm[6])
                    else:
                    #TODO fix the injection
                    #unstable, turned off in Main core
                        outputs = FastEuler.euler_m_slave( sldriver.Ib, sldriver.Ip,slfiltered,sldriver.Dt,sldriver.L,sldriver.nPointsPerPulse,sldriver.time,output[1],self.widgets.injfix.var.get(),slavearg,cm[5],cm[4],cm[7],0,0,cm[2],cm[3],cm[6])
                        
                        #noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None):   
                        #cm[5],cm[4],cm[7],0,0,cm[2],cm[3],cm[6])                        
                        
                        outputs[1] = np.sqrt(outputs[1]) * np.exp(1j*outputs[3])
                if(self.widgets.option_var.get().find('WDM') != -1 ):
                    outputs[1] = wdm.wdm_sigfilter( wdm, outputs[1], float(self.slavewdmdata[0]), int(self.slavewdmdata[1]), float(self.slavewdmdata[2]),str(self.widgets.wdmvar2.get()))
            
            
            
            if(self.slavestate == 1):
                detector_outS = wdm.photodetector( wdm, outputs[1] , float(self.photdata[0]), int(self.photdata[1]), float(self.photdata[2]),"Bessel")
                if(self.widgets.power_output.var.get()):   
                    
                    Qprev = detector_outS
                    #hbar = 1.054 *1e-25 #J*s
                    #t = 0.001 #*1e-9 #s                 
                    #lambd = 1550#*1e-9 #m #1550 #nm
                    c = 2*1e17# m/s
                    #omega = c / laserS.lambd*1e3
                    detector_outS = ( laserS.eta*laserS.hbar *(c/(laserS.lambd * 1e3) ) )/(2*laserS.Gamma*laserS.tau_ph)* Qprev                     

            phiout = np.zeros(len(output[0]))
            nuout = np.zeros(len(output[0]))
            interferout = np.zeros(len(output[0]))
              
            if (cm[2] == 1):
                phiout = output[3]                      #getting phase output in the plot array
                
            if (cm[3] == 1 ):
                if ( self.widgets.scipym.var.get() and self.widgets.calc_var.get() == "Runge_Kutta"):
                    nuout = nu#output[4]
                else:
                    nuout= output[4]
                
                #TODO!
            if (cm[7] == 1):
                #interferout = laser.interference(output[1])           
                interferout = laser.interference(output[1])#( np.sqrt(output[1])  * np.exp(1j*output[3])  )
                
            if (dm[0] == 1):                            #adding calculation results to the plot array
                outputarr.append(filtered*1e12)
                outputlabels.append(r"Master Electric Pulses, $mA$")
            if (dm[1] == 1):                            
                outputarr.append(output[0])
                outputlabels.append("Master Carrier Number")
            if (dm[2] == 1):
                outputarr.append(detector_outM)
                #outputarr.append(output[1]) 
                if (self.widgets.power_output.var.get()):
                    outputlabels.append("Master Intensity, $W$") 
                else:
                    outputlabels.append("Master Intensity")
            if (dm[3] == 1):                
                outputarr.append(phiout)
                outputlabels.append("Master Phase")
            if (dm[4] == 1):
                outputarr.append(nuout)
                outputlabels.append(r"Master Frequency, $Ghz$")
            if (dm[5] == 1):
                 outputarr.append(output[2])
                 outputlabels.append(r"Master Temperature variation, $K$")
            if (dm[6] == 1):
                outputarr.append(interferout)
                outputlabels.append("Master Interference")
                
                
                
            if (self.widgets.spectrumreg_on.var.get()):
                text_spec = []
                self.tmp_spec = []
                if (self.widgets.spectrum_menu_reg.option_var.get() == "Master"):
                    arrsize = int(2*len(output[1])-1)
                    #print(f"arrsize {arrsize} {driver.nPointsPerPulse} {len(output[1])}")
                    samplingFrequency   = 1 / driver.Dt
                    freq = np.linspace(-1/2, 1/2,arrsize)*samplingFrequency
                    norm_pulses =np.zeros(arrsize, dtype = complex)
                    for i in range(0,int(driver.nPointsPerPulse)):
                        norm_pulses[i] = output[1][i]
                
                    Spectrum = abs(np.fft.fftshift(np.fft.fft(norm_pulses) ))**2
                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))
                    axs.set_title("Spectrum")
                    axs.set_xlim([-1000, 1000])
                    axs.set_yscale('log')                
                    axs.grid(which = 'both' , linewidth =1)
                    axs.set(xlabel=r'$Ghz$', ylabel="")
                    axs.plot(freq,Spectrum)
                    fig.tight_layout()
                    plt.show() 
                    text_spec = np.c_[np.stack([freq, Spectrum], axis = 1)]
                else:
                    if(self.slavestate == 1): 
                        arrsize = int(2*len(outputs[1])-1)
                        samplingFrequency   = 1 / driver.Dt
                        freq = np.linspace(-1/2, 1/2,arrsize)*samplingFrequency
                        norm_pulses =np.zeros(arrsize, dtype = complex)
                        for i in range(0,int(sldriver.nPointsPerPulse)):
                            norm_pulses[i] = outputs[1][i]
                    
                        Spectrum = abs(np.fft.fftshift( np.fft.fft( norm_pulses)  ) )**2  
                        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))
                        axs.set_title("Spectrum")
                        axs.set_xlim([-1000, 1000])
                        axs.set_yscale('log')                
                        axs.grid(which = 'both' , linewidth =1)
                        axs.set(xlabel=r'$Ghz$', ylabel="")
                        axs.plot(freq,Spectrum)
                        fig.tight_layout()
                        plt.show() 
                        text_spec = np.c_[np.stack([freq, Spectrum], axis = 1)]
                
                self.tmp_spec = text_spec
                np.savetxt('results/last_spectrum_calculation_regular.txt', text_spec)
                self.widgets.download_spec.place(x = 250, y= 750, width = 200, height = 50)
            else:
                self.widgets.download_spec.place_forget()
            if(self.slavestate == 1):                   #adding slave calculation results to the plot array
                phiSout = np.zeros(len(outputs[0]))
                nuSout = np.zeros(len(outputs[0]))
                interferSout = np.zeros(len(outputs[0]))
                # if (cs[2] == 1):
                #     phiSout = outputs[3]
      
                if (cs[3] == 1):
                    if ( self.widgets.scipym.var.get() and self.widgets.calc_var.get() == "Runge_Kutta"):
                        nuSout = nuS#output[4]
                    else:
                        nuSout = outputs[4]
                if (cs[7] == 1):
                    #interferSout = FastEuler.interference(outputs[1], laserS.L, laserS.Dt, laserS.delay, laserS.phase_m1, laserS.phase_m2)
                    interferSout = laserS.interference(outputs[1])
                print("Slave is active")
                if( ds[0] == 1):
                    outputarr.append(slfiltered*1e12)
                    outputlabels.append(r"Slave Electric Pulses, $mA$")
                if( ds[1] == 1):
                    outputarr.append(outputs[0])
                    outputlabels.append("Slave Carrier Number")
                if( ds[2] == 1):
                    #outputarr.append(outputs[0])
                    outputarr.append(detector_outS)
                    #outputlabels.append("Slave Carrier Number")
                    if(self.widgets.power_output.var.get()):
                        outputlabels.append(r"Slave Intensity, $W$")
                    else:
                        outputlabels.append("Slave Intensity")
                if( ds[3] == 1):
                    outputarr.append(outputs[3])
                    outputlabels.append("Slave Phase")
                if (ds[4] == 1):
                    outputarr.append(nuSout)
                    outputlabels.append("Slave Frequency, $Ghz$")
                if (ds[5] == 1):
                    if ( self.widgets.scipym.var.get()):
                        #nuSout = nuS#output[4]
                        outputarr.append(outputs[2])
                    else:
                        outputarr.append(outputs[2])
                    outputlabels.append(r"Slave Temperature variation, $K$")
                if (ds[6] == 1):
                    outputarr.append(interferSout)
                    outputlabels.append("Slave Interference")
                 
            if (cm[8] == 1):      #stability diagram, displayed in a separate window 
                if (mIb > laser.Ith):
                    
                    x = np.linspace(float(self.transdata[0]), float(self.transdata[1]))
                    plt.yscale("log")
                    plt.xscale("log")
                    plt.xlabel(r"$\omega_p$")
                    plt.ylabel(r"$H(i*w_p)$")
                    responsedata = laser.transfer(x,mIb)
                    tmax = np.ndarray.max(responsedata)
                    tmax_i = np.ndarray.argmax(responsedata)
                    tmax = round(tmax,4)
                    
                    #print(str(tmax) + " " + str(responsedata[tmax_i]) + " maximum value on stability diagram")
                    #self.widgets.stabbutton.config(command = Driver.getfig_n(self,1,output[0],output[1],"","Detuning frequency",self.parent, self.stabcanvas, self.stabfigure,0,self.widgets.pulseframe_stab))
                    
                    plt.plot(x, responsedata, color = 'red',markersize = 3)
                    plt.plot(x[tmax_i],responsedata[tmax_i], '^')
                    plt.text(x[tmax_i], responsedata[tmax_i]*2, "maximum at : $\omega_\pi$= " + str(round(x[tmax_i],2))+ "  H= " + str(round(responsedata[tmax_i],2)) , style='italic',
                            bbox={'facecolor': 'red', 'alpha': 0.8, 'pad': 8})
                    #plt.grid(which = 'both' , linewidth =1)
                    plt.yticks([0.1,0.25,1,10])
                    plt.xticks([1,float(self.transdata[1])/10,float(self.transdata[1] )])
                    
                    plt.show()
                    tr_out =[]
                    tr_out.append(x)                    
                    tr_out.append(responsedata)
                    self.tnp = np.c_[np.stack(tr_out, axis = 1)]  
                    np.savetxt('results/last_transfer.txt', self.tnp)
                else:
                    print(" Response function : Current value is too low")

  
            text_out =[]
            text_out.append(maxtime)
            for i in range(0,len(outputarr)):               
                text_out.append(outputarr[i])
            self.tnp = np.c_[np.stack(text_out, axis = 1)]  
            np.savetxt('results/last_calculation.txt', self.tnp)
             
            for i in range(0,len(outputarr)):
                outputleg.append("")


            if (self.widgets.graphic_output.var.get()):
                color = ["tab:orange","tab:green","tab:red","tab:blue"]
                left  = 0.118  # the left side of the subplots of the figure
                right = 0.328    # the right side of the subplots of the figure
                bottom = 0.764   # the bottom of the subplots of the figure
                top =0.94      # the top of the subplots of the figure
                wspace = 0   # the amount of width reserved for blank space between subplots
                hspace = 0.425   # the amount of height reserved for white space between subplots
            
                fig, axs = plt.subplots(nrows=len(outputarr), ncols=1, figsize=(16, 9))
                #fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 3))
                fig.subplots_adjust(left =left , bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
                print(f"len {len(maxtime)} + lenout {len(outputarr[0])} + {len(name)} ")
                if(len(outputarr) > 1):
                    for i in range(0,len(outputarr)): 
                    
                        if (outputlabels[i].find('Interference') != -1 ):
                            shift = int(laser.delay / laser.Dt)
                            axs[i].plot(maxtime[0:len(outputarr[i])-shift],outputarr[i][0:len(outputarr[i])-shift],color[i%len(color)])
                        else:
                            axs[i].plot(maxtime,outputarr[i],color[i%len(color)])
                       # print(str(name[i]) + " " + str(len(y[i])) +" " +  str(len(x[0:len(y[i])])))
                        axs[i].set_title(outputlabels[i]) 
                        axs[i].yaxis.major.formatter._useMathText = True
                        if ( i == len(outputarr)-1):
                            axs[i].set(xlabel='T ns', ylabel="")
                else:                                     
                     axs.plot(maxtime,outputarr[0])#,"go",markersize = 3)
                     axs.set_title(outputlabels[0])
                     axs.yaxis.major.formatter._useMathText = True
                     axs.set(xlabel='T ns', ylabel="")
                     axs.grid(which = 'both' , linewidth =1)
                fig.tight_layout()
                plt.show()
            
            
            toolbar = NavigationToolbar2Tk(self.canvas, self.widgets.pulseframe_reg, pack_toolbar=False)
            toolbar.update()
            self.parent.update()
            
            inter_delay = 0
            if(dm[6] == 1 or ds[6] == 1):
                inter_delay = int(laser.delay / laser.Dt)
                print(f'inter delay : {inter_delay}')
            # if ( self.widgets.scipym.var.get()):
            #     self.widgets.button.config(command = Driver.getfig_n(self,toolbar, len(outputarr), maxtime,outputarr , outputlabels, outputleg, self.parent, self.canvas, self.figure,shift = inter_delay, frame = self.widgets.pulseframe_reg,W = 650,H = 600))  
            # else:
            self.widgets.button.config(command = Driver.getfig_n(self,toolbar, len(outputarr), maxtime,outputarr , outputlabels, outputleg, self.parent, self.canvas, self.figure,shift = inter_delay, frame = self.widgets.pulseframe_reg, W = 650, H =600))  
            print("Time" , datetime.now() - start)
        #self.parent.update()
    def arbpulsecalculation(self): #perform calculation in the "Modified Pulses" frame
        
        cm = [self.widgets.m_checkbox_clc[i].var.get() for i in range(0,self.widgets.master_clcoptions)]
        dm = [self.widgets.m_checkbox_plt[i].var.get() for i in range(0,self.widgets.master_pltoptions)]
        #cm - calculation master option, dm - display slave option
        #self.modfigure = plt.Figure(figsize=(16,9), dpi=100)    #plot display for arbitrary pulses
        #self.modcanvas = FigureCanvasTkAgg(self.modfigure, self.widgets.pulseframe_mod)
        self.savelaser()
        mIb = float(self.data[0])*1e-12
        mIp = float(self.data[1])*1e-12
        mrate = float(self.data[3])
        mbandwidth = float(self.data[4])
        mforder = int(self.data[5])
        mfp = float(self.data[6])
        mw = float(self.data[7])        
        mclctime = float(self.data[8])
        
        arbfinaloutput = []
        sarbfinaloutput = []
        arbfinaltime = []
        arbleg = []
        arbname =[]
        arbplotdata = [] #stores data going to the plot module 
        if(self.pulseid == 0):
            print("No arbitrary pulses from master laser")
            return 0
            
        if ( self.data[2] != ""):  
            npulses = float(self.data[2] )
        mrate = self.Arbpulses[0].data[8]
        for i in range(0,self.pulseid):                         #creating arbitraray pulses from main lasers
            for j in range(0,10):
                if (self.Arbpulses[i].entry[j].get() != ''):
                    self.Arbpulses[i].data[j] = self.Arbpulses[i].entry[j].get() 
        if(self.pulseid > 0):                                   #arbitrary pulses share repetition rate and lc time
            for i in range(0,self.pulseid):
                self.Arbpulses[i].data[6] = self.Arbpulses[0].data[6] #REP RATE
                self.Arbpulses[i].data[8] = self.Arbpulses[0].data[8] #INTEGRATION RATE
                self.Arbpulses[i].data[7] = self.Arbpulses[0].data[7] #TIME SHIFT
                #self.Arbpulses[i].data[9] = self.Arbpulses[0].data[9] #INCLINE
            totalfp = 0
            totalw = 0
            for i in range(0,self.pulseid):                     #calculate total repetition rate and width across all pulses
                totalfp = totalfp + 1/self.Arbpulses[i].data[6]
                totalw = totalw + self.Arbpulses[i].data[2]
            totalfp = round(totalfp,4)
            totalw =  round(totalw,4)
            print("total fp " + str(totalfp) + "totalw " + str(totalw))

            
            Dt = self.Arbpulses[0].data[8]
            time = np.arange(0,float(totalfp),Dt )
            print("Total time for all master pulses " +str(len(time)))
            
            maxfp   = totalfp
            maxtime = time
            start = datetime.now() 
            driverM =  Driver(self.Arbpulses[i].data[1]*1e-12,self.Arbpulses[i].data[0]*1e-12,self.Arbpulses[0].data[6],self.Arbpulses[i].data[2],mclctime,self.pulseid,self.Arbpulses[0].data[8])
            print("MAIN DRIVER " +str(len(driverM.time)) )
          
      
            sigtype = self.widgets.sigtype_menu.option_var.get()
            for i in range(0,self.pulseid):
                #print(f'current incline is {self.Arbpulses[i].data[9]}')
                self.Arbpulses[i].data[6] = self.Arbpulses[0].data[6]
                driver =  Driver(mIp,mIb,self.Arbpulses[i].data[6],self.Arbpulses[i].data[2],mclctime,self.pulseid,self.Arbpulses[0].data[8])
                arboutput = driver.arbpulses(driver,self.Arbpulses[i].data[0]*1e-12,self.Arbpulses[i].data[1]*1e-12,self.Arbpulses[i].data[2],\
                self.Arbpulses[i].data[3]*1e-12,self.Arbpulses[i].data[4],self.Arbpulses[i].data[5],self.Arbpulses[i].data[6],self.Arbpulses[i].data[7],sigtype,self.Arbpulses[i].data[9])
                                                
                arbmin= driver.arbminor(self,self.Arbpulses[i].data[2],self.Arbpulses[i].data[3]*1e-12,self.Arbpulses[i].data[4],self.Arbpulses[i].data[5],self.Arbpulses[i].data[6],0,0,"Bessel")
                arbmin = arbmin[0:len(arboutput)]
                arboutput = arboutput - arbmin[0:len(arboutput)]
                
                arbfinaloutput = np.concatenate([arbfinaloutput,arboutput])
                #arbfinaloutput = np.concatenate( [arbfinaloutput , arbfiltered])
            arbfiltered = driver.arbbetterfilter(self.Arbpulses[i].data[0]*1e-12, self.Arbpulses[i].data[1]*1e-12, self.Arbpulses[i].data[3]*1e-12, self.Arbpulses[i].data[4], self.Arbpulses[i].data[5], arbfinaloutput, mbandwidth, mforder, "Bessel")

            if (self.slavepulseid > 0):
                driverS =  Driver(self.slaveArbpulses[i].data[1]*1e-12,self.slaveArbpulses[i].data[0]*1e-12,self.Arbpulses[0].data[6],self.slaveArbpulses[i].data[2],mclctime,self.slavepulseid,self.Arbpulses[0].data[8])
                print("SLAVE DRIVER " +str(len(driverS.time)) )
                for i in range(0,self.slavepulseid):
                    for j in range(0,10):
                        if (self.slaveArbpulses[i].entry[j].get() != ''):
                            self.slaveArbpulses[i].data[j] = self.slaveArbpulses[i].entry[j].get() 
                    
                for i in range(0,self.slavepulseid):
                    self.slaveArbpulses[i].data[6] = self.slaveArbpulses[0].data[6]
                    self.slaveArbpulses[i].data[7] = self.slaveArbpulses[0].data[7]
                    self.slaveArbpulses[i].data[9] = self.slaveArbpulses[0].data[9]
                stotalfp = 0
                stotalw = 0
                for i in range(0,self.slavepulseid):
                    stotalfp = stotalfp + 1/self.slaveArbpulses[i].data[6]
                    stotalw = stotalw + self.slaveArbpulses[i].data[2]
                stotalfp = round(stotalfp,4)
                stotalw =  round(stotalw,4)
                print("slave total fp " + str(stotalfp) + "slave totalw " + str(stotalw))

                Dt = self.Arbpulses[0].data[8]
                slavetime = np.arange(0,float(stotalfp),Dt ) 
                print("Total time for all slave pulses " +str(len(time)))
                for i in range(0,self.slavepulseid):
                    self.slaveArbpulses[i].data[6] = self.slaveArbpulses[0].data[6]
                    driver =  Driver(mIp,mIb,self.slaveArbpulses[i].data[6],self.slaveArbpulses[i].data[2],mclctime,self.pulseid,self.Arbpulses[0].data[8])
                    arboutput = driver.arbpulses(driver,self.slaveArbpulses[i].data[0]*1e-12,self.slaveArbpulses[i].data[1]*1e-12,self.slaveArbpulses[i].data[2],\
                    self.slaveArbpulses[i].data[3]*1e-12,self.slaveArbpulses[i].data[4],self.slaveArbpulses[i].data[5],self.slaveArbpulses[i].data[6],self.slaveArbpulses[i].data[7],sigtype,self.slaveArbpulses[i].data[9])
                    
                    arbmin= driver.arbminor(self,self.slaveArbpulses[i].data[2],self.slaveArbpulses[i].data[3]*1e-12,self.slaveArbpulses[i].data[4],self.slaveArbpulses[i].data[5],self.slaveArbpulses[i].data[6],0,0,"Bessel")
                    arbmin = arbmin[0:len(arboutput)]
                    arboutput = arboutput - arbmin[0:len(arboutput)]
                    sarbfinaloutput =np.concatenate([sarbfinaloutput,arboutput])
                sarbfiltered = driver.arbbetterfilter(self.slaveArbpulses[i].data[0]*1e-12, self.slaveArbpulses[i].data[1]*1e-12, self.slaveArbpulses[i].data[3]*1e-12, self.slaveArbpulses[i].data[4], self.slaveArbpulses[i].data[5], sarbfinaloutput, mbandwidth, mforder, "Bessel") 

                if(totalfp < stotalfp):
                    maxtime = slavetime
                    maxfp = stotalfp
                elif(totalfp > stotalfp):
                    maxtime = time
                    maxfp = totalfp
                if (self.pulseid < self.slavepulseid ): 
                    arbfiltered =np.concatenate( [arbfiltered, np.zeros(abs(len(arbfiltered) - len(maxtime))) ])
                    driverM.time = driverS.time#np.concatenate( [driverM.time, np.zeros(abs(len(driverM.time) - len(maxtime))) ],axis = 0)
                    driverM.L = len(maxtime)
                elif (self.pulseid > self.slavepulseid ): 
                    sarbfiltered = np.concatenate([sarbfiltered , np.zeros(abs(len(sarbfiltered) - len(maxtime))) ])
                    driverS.time = driverM.time#np.concatenate( [driverS.time, np.zeros(abs(len(driverS.time) - len(maxtime))) ],axis = 0)
                    driverS.L = len(maxtime)
                
                #startshift = self.slaveArbpulses[0].data[7]
               #if (startshift > 0):
                #    startshiftadd = np.zeros(int(startshift/Dt))
                 #   sarbfiltered=np.concatenate([startshiftadd,sarbfiltered])
                  #  sarbfiltered = sarbfiltered[0:len(maxtime)]
                #slaveupperTime = self.pulseid/fp
            print("Highest total fp " + str(maxfp) +" " + str(len(maxtime)) )
            
            #electric pulse preparation is done
            arg = []
            slavearg = []
            for i in range(self.datalstart,self.datalend):
                arg.append(float(self.data[i]))
            for i in range(self.datalstart,self.datalend):
                slavearg.append(float(self.slavedata[i]))
            arg[25] = float(self.intdata[0])   #interferometer delay
            arg[26] = float(self.intdata[1])
            arg[27] = float(self.intdata[2])
            slavearg[25] = arg[25]
            slavearg[26] = arg[26]
            slavearg[27] = arg[27]
                 
            laser = Laser(mIp,driverM,arg)
            #print(str(len(filtered)) + " " + str(len(slfiltered)) +  " "+ str(len(driver.time)) + " " + str(len(sldriver.time)))
            outputarr = []
            outputlabels = []
            outputleg = []
            timearr = []

            name = ""
            print(len(arbfiltered))
            
            
            if ( self.widgets.calc_var.get() == "Runge_Kutta"):     
                
                if ( self.widgets.scipym.var.get()):
                    output = laser.sci_ODE(arbfiltered,cm[2],cm[3],cm[5],cm[4],cm[7])
                    print(f"inj scipy {len(output[1])}")
                    tempphase = output[2]
                    tempT = output[3]
                    # output[2] = tempT
                    # output[3] = tempphase
                    nu = np.zeros(len(output[4]))
                    for i in range(0,len(output[4])-1):
                        nu[i] = (output[3][i+1] - output[3][i]) /(laser.Dt*2*np.pi)
                        
                else:
                    output  = laser.runge_kutta(arbfiltered,cm[2],cm[3],cm[5],cm[4],cm[7])#laser.euler_m(arbfiltered,cm[2],cm[3],cm[5],cm[4],cm[7])
            else:                
                
                #self,I,phase = None,freq = None, noise = None, temp = None, interference = None, field = None, biexp = None
                if (self.widgets.Euler_mod.option_var.get() == 'Regular'):
                    output = laser.euler_m(arbfiltered,cm[2] , cm[3], cm[5] , cm[4],cm[7],0,0)
                else:
                    output = FastEuler.fast_euler( driverM.Ib, driverM.Ip,arbfiltered,driverM.Dt,driverM.L,driverM.nPointsPerPulse,arg,cm[2],cm[3],cm[5],cm[4],cm[7],0,0)
                    output[1] = np.sqrt(output[1]) * np.exp(1j*output[3])
                    
                #output= FastEuler.fast_euler( driverM.Ib, driverM.Ip,arbfiltered,driverM.Dt,driverM.L,driverM.nPointsPerPulse,arg,cm[2],cm[3],cm[5],cm[4],cm[7],0,0)
                
                # c0 = self.widgets.m_checkbox_clc[0].var.get()    #driver
                # c1 = self.widgets.m_checkbox_clc[1].var.get()    #laser output  N Q
                # c2 = self.widgets.m_checkbox_clc[2].var.get()    #phase
                # c3 = self.widgets.m_checkbox_clc[3].var.get()    #frequency 
                # c4 = self.widgets.m_checkbox_clc[4].var.get()    #temperature
                # c5 = self.widgets.m_checkbox_clc[5].var.get()    #noise
                # c6 = self.widgets.m_checkbox_clc[6].var.get()    #biexp temp
                # c7 = self.widgets.m_checkbox_clc[7].var.get()    #interference
                
                #0.91 version
                # c0 = self.cdata[0]    #driver
                # c1 = self.cdata[1]    #laser output  N Q
                # c2 = self.cdata[2]    #phase
                # c3 = self.cdata[3]    #frequency 
                # c4 = self.cdata[4]    #complex field
                # c5 = self.cdata[5]    #temperature
                # c6 = self.cdata[6]    #noise
                # c7 = self.cdata[7]    #biexp temp
                # c8 = self.cdata[8]    #interference
                #output  = laser.euler_m(laser,arbfiltered,cm[5],cm[4],cm[7],0,0,cm[2],cm[3],0)
                
                
                
           #print ("MASTER INJECTION PARAMETERS: kappa_inj " + str(laser.kap_inj) + " w_inj " + str(laser.w_inj))
            #print ("MASTER INJECTION PARAMETERS: w_M " + str(laser.w_M) + " om_th " + str(laser.om_th))
            wdm = Filter(mfp,mw,mclctime,npulses,self.Arbpulses[0].data[8]) 
            if(self.widgets.option_var.get().find('Master-WDM') != -1 ): 
                #print(self.wdmdata[0] + " " + self.wdmdata[1] + " " + self.wdmdata[2] + " wdm data " + str(self.widgets.wdmvar1.get()))
                output[1] = wdm.wdm_sigfilter( wdm, output[1], float(self.wdmdata[0]), int(self.wdmdata[1]), float(self.wdmdata[2]),str(self.widgets.wdmvar1.get()))
           
            detector_outM = wdm.photodetector( wdm, output[1] , float(self.photdata[0]), int(self.photdata[1]), float(self.photdata[2]),"Bessel")
            #plt.plot(maxtime, output[0],linewidth=2,color = "tab:orange")
            

            ##################################### slave driver
            slIb = float(self.slavedata[0])*1e-12
            slIp = float(self.slavedata[1])*1e-12


            cs = [self.widgets.s_checkbox_clc[i].var.get() for i in range(0,self.widgets.slave_clcoptions)]
            ds  = [self.widgets.s_checkbox_plt[i].var.get() for i in range(0,self.widgets.slave_pltoptions)]      
            #cs - calculation slave options, ds - display slave options
                        

            if(self.slavestate == 1 and self.slavepulseid >0):
                laserS = Laser(mIp,driverS,slavearg)
                #self,EM,driver,I,noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None):  
                #outputs = laserS.euler_m_slave(output[1],driverS,sarbfiltered,cs[5],cs[4],cs[7],0,0,cs[2],cs[3],0)
                if ( self.widgets.calc_var.get() == "Runge_Kutta"):              
                    
                    if ( self.widgets.scipym.var.get()):   
                       #unstable, turned off
                       outputs = laserS.sci_ODE_slave(sarbfiltered,output[1],cs[2],cs[3],cs[5],cs[4],cs[7],cs[6])
                       nuS = np.zeros(len(outputs[4]))
                       for i in range(0,len(outputs[3])-1):
                           nuS[i] = (outputs[3][i+1] - outputs[3][i]) /(laser.Dt*2*np.pi)
                    else:       
                        outputs = laserS.runge_kutta_slave(output[1], sarbfiltered,self.widgets.injfix.var.get(),cs[2],cs[3],cm[5],cm[4],cm[7],cm[6])
                    #def runge_kutta_slave(self,EM,I,injfix, phase = None, freq =None, noise = None, temp = None, interference = None, field = None):
                else:
                    #outputs = FastEuler.euler_m_slave( driverS.Ib, driverS.Ip,sarbfiltered,driverS.Dt,driverS.L,driverS.nPointsPerPulse,driverS.time,output[1],self.widgets.injfix.var.get(),slavearg,cm[2],cm[3],cm[5],cm[4],cm[7],0,0)
                    
                    #actual
                    #outputs = laserS.euler_m_slave(output[1],driverS,sarbfiltered,self.widgets.injfix.var.get(),cs[5],cs[4],cs[7],0,0,cs[2],cs[3],cm[6])
  
    #TODO!              outputs = laserS.euler_m_slave(output[1],sldriver,slfiltered,self.widgets.injfix.var.get(),cs[5],cs[4],cs[7],0,0,cs[2],cs[3],0)
                    if (self.widgets.Euler_mod.option_var.get() == 'Regular'):
                         outputs = laserS.euler_m_slave(output[1],driverS,sarbfiltered,self.widgets.injfix.var.get(),cs[5],cs[4],cs[7],0,0,cs[2],cs[3],cm[6])
                    else:
                    #TODO fix the  injection
                        #print(f"{driverS.Ib}, {driverS.Ip},{driverS.Dt},{driverS.L}, {driverS.nPointsPerPulse}")
                        #print(output[1])
                        #unstable, turned off
                        outputs = FastEuler.euler_m_slave( driverS.Ib, driverS.Ip,sarbfiltered,driverS.Dt,driverS.L, driverS.nPointsPerPulse, driverS.time,output[1],self.widgets.injfix.var.get(),slavearg,cm[5],cm[4],cm[7],0,0,cm[2],cm[3],cm[6])
                        
                        #noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None):   
                        #cm[5],cm[4],cm[7],0,0,cm[2],cm[3],cm[6])

                        outputs[1] = np.sqrt(outputs[1]) * np.exp(1j*outputs[3])  
  
    
                  #0.91
                #   for i in range(17,26): #same calc parameters for slave laser
                # cs[i] = self.cdata[i]
                #self,EM,driver,I,noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None
                
                # c0 = self.widgets.m_checkbox_clc[0].var.get()    #driver
                # c1 = self.widgets.m_checkbox_clc[1].var.get()    #laser output  N Q
                # c2 = self.widgets.m_checkbox_clc[2].var.get()    #phase
                # c3 = self.widgets.m_checkbox_clc[3].var.get()   #frequency 
                # c4 = self.widgets.m_checkbox_clc[4].var.get()    #temperature
                # c5 = self.widgets.m_checkbox_clc[5].var.get()    #noise
                # c6 = self.widgets.m_checkbox_clc[6].var.get()    #biexp temp
                # c7 = self.widgets.m_checkbox_clc[7].var.get()    #interference
                
                  #outputs = laserS.euler_m_slave(output[1],driverS,sarbfiltered,cs[5],cs[4],cs[7],0,0,cs[2],cs[3],0)
                  
                  
                  
                if(self.widgets.option_var.get().find('Slave-WDM') != -1 ):
                    print(self.slavewdmdata[0] + " " + self.slavewdmdata[1] + " " + self.slavewdmdata[2] + " slave wdm data " + str(self.widgets.wdmvar1.get()))
                    outputs[1] = wdm.wdm_sigfilter( wdm, outputs[1], float(self.slavewdmdata[0]), int(self.slavewdmdata[1]), float(self.slavewdmdata[2]),str(self.widgets.wdmvar2.get()))
                
                #if(self.widgets.option_var.get().find('Slave-WDM') != -1 ):
                #    print(self.slavewdmdata[0] + " " + self.slavewdmdata[1] + " " + self.slavewdmdata[2] + " slave wdm data " + str(self.widgets.wdmvar1.get()))
                #    outputs[1] = wdm.wdm_sigfilter( wdm, outputs[1], float(self.slavewdmdata[0]), int(self.slavewdmdata[1]), float(self.slavewdmdata[2]),str(self.widgets.wdmvar2.get()))
                
                #if(self.slavestate == 1):
                detector_outS = wdm.photodetector( wdm, outputs[1] , float(self.photdata[0]), int(self.photdata[1]), float(self.photdata[2]),"Bessel")
                
                if(self.widgets.power_output.var.get()):   
                    
                
                    c = 2*1e17# m/s
                    omega = c / laserS.lambd* 1e3
                    # P = eta  * hbar * omega * Q / (2 * G * t)    
                    #detector_outM = ( laser.eta*hbar *(c/laser.lambd * 1e3 ) )/(2*laser.Gamma*laser.tau_ph)* Q                    
                    Qprev = detector_outS
                    detector_outS = laserS.eta*laserS.hbar *(c/(laserS.lambd * 1e3) )/(2*laserS.Gamma*laserS.tau_ph)*Qprev
                    
            if (self.widgets.power_output.var.get()):
                Q = detector_outM                
                c = 2*1e17# m/s                
                detector_outM = ( laser.eta*laser.hbar *(c/(laser.lambd * 1e3) ) )/(2*laser.Gamma*laser.tau_ph)* Q
                
            phiout = np.zeros(len(output[0]))
            nuout = np.zeros(len(output[0]))
            #interferout = np.zeros(len(output[0]))
            
            if (cm[2] == 1):
                phiout = output[3]                      #getting phase output in the plot array
                
                
            if (cm[3] == 1):
                if ( self.widgets.scipym.var.get() and self.widgets.calc_var.get() == "Runge_Kutta"):
                    nuout = nu#output[4]
                else:
                    nuout= output[4]
                    

            if (cm[7] == 1):
                interferout = laser.interference(output[1])
            #print( "slave size " + str(len(outputs[0])))
            
            phiSout = np.zeros(len(output[0]))
            nuSout = np.zeros(len(output[0]))
            interferSout = np.zeros(len(output[0]))
           
            if (dm[0] == 1):
                arbleg.append(" ")
                arbname.append(r"Master Electric Pulses, $mA$")
                arbplotdata.append(arbfiltered*1e12)
            if (dm[1] == 1):
                

                arbname.append("Master Carrier number") 
                arbleg.append(" ")
                arbplotdata.append(output[0])
                    
            if (dm[2] == 1):
                
                if (self.widgets.power_output.var.get()):
                    arbname.append("Master Intensity, $W$") 
                else:
                    arbname.append("Master Intensity")
                arbleg.append(" ")
                #arbplotdata.append(output[0])
                arbplotdata.append(detector_outM)
            if (dm[3] == 1):                
                #phiout = laser.phasepack(output[1])
                arbplotdata.append(phiout)
                arbname.append("Master Phase ")
                arbleg.append(" ")
                #text_output.append(phiout)
            if (dm[4] == 1):
                #nuout = laser.frequency(phiout)
                arbplotdata.append(nuout)
                arbname.append("Master Frequency, $Ghz$")
                arbleg.append(" ")
            if (dm[5] == 1):
                arbplotdata.append(output[2])
                arbname.append(r"Master Temperature difference, $K$")
                arbleg.append(" ")    
            if (dm[6] == 1):
                #interferout = laser.interference(output[1])
                arbplotdata.append(laser.interference(output[1]))
                arbleg.append(" ")
                arbname.append("Master Interference")
                
            if(self.slavestate == 1 and self.slavepulseid > 0):   
                phiSout = np.zeros(len(outputs[0]))
                nuSout = np.zeros(len(outputs[0]))
                interferSout = np.zeros(len(outputs[0]))
                if (cs[2] == 1):
                    phiSout = outputs[3]#laserS.phasepack(outputs[1])
                # if (cs[3] == 1):
                #     nuSout = outputs[4]#laserS.frequency(phiSout)
                    
                if (cs[3] == 1):
                    if ( self.widgets.scipym.var.get() and self.widgets.calc_var.get() == "Runge_Kutta"):
                        nuSout = nuS#output[4]
                    else:
                        nuSout = outputs[4]
                        
                if (cs[7] == 1):
                    interferout = laserS.interference(outputs[1])
                print("Slave is active")
                    
                if(ds[0] == 1):
                    arbleg.append(" ")
                    arbname.append(r"Slave Electric pulses, $mA$")
                    arbplotdata.append(sarbfiltered*1e12)
                if (ds[1] == 1):
                    arbname.append("Slave Carrier number") 
                    arbleg.append(" ")
                    arbplotdata.append(outputs[0])
                
                if( ds[2] == 1):
                    
                    if(self.widgets.power_output.var.get()):
                        arbname.append(r"Slave Intensity, $W$")
                    else:
                        arbname.append("Slave Intensity")
                        
                    #arbname.append("Slave Intensity")
                    #outputarr.append(outputs[0])
                    arbplotdata.append(detector_outS)
                    arbleg.append(" ")
                if( ds[3] == 1):
                    
                    #phiSout = laserS.phasepack(outputs[1])
                    arbplotdata.append(phiSout)
                    arbleg.append(" ")
                    arbname.append("Slave Phase")
                if (ds[4] == 1):

                    arbplotdata.append(nuSout)
                    arbleg.append(" ")
                    arbname.append("Slave Frequency, $Ghz$")
                    #text_output.append(outputs[3])
                if (ds[5] == 1):             
                    arbplotdata.append(outputs[2])
                    arbname.append(r"Slave Temperature Difference, $K$")
                    arbleg.append(" ")
                if (ds[6] == 1):
                    interferSout = laserS.interference(outputs[1])
                    arbplotdata.append(interferSout)
                    arbname.append("Slave Interference")
                    arbleg.append(" ")
                    

            text_out =[]
            text_out.append(maxtime)
            for i in range(0,len(arbplotdata)):
                text_out.append(arbplotdata[i])
            self.tnp = np.c_[np.stack(text_out, axis = 1)]  
            np.savetxt('results/last_arbpulse_calculation.txt', self.tnp)
            
            self.widgets.threscurrM.var.set(round(laser.Ith*1e12,4))
            if(self.slavestate == 1 and self.slavepulseid >0):
                self.widgets.threscurrS.var.set(round(laserS.Ith*1e12,4))

            
            
            #arbsignal = driver.arbpulses(driver,self.arbdata[0]*1e-12,self.arbdata[1]*1e-12,self.arbdata[2],self.arbdata[3]*1e-12,self.arbdata[4],self.arbdata[5],self.arbdata[6])
            #arbfiltered = driver.arbfilter(driver, self.arbdata[0]*1e-12, self.arbdata[1]*1e-12, arbsignal, 10, 2, 'Bessel')
            
            pulse1 = self.widgets.specpulse1.var.get()
            pulse2 = self.widgets.specpulse2.var.get()
            
            if (self.widgets.spectrumall_on.var.get()):
                #print(f" S driver {driverS.Dt} {driverS.nPointsPerPulse} {len(driverS.time)}")
                #print(f" M driver {driver.Dt} {driver.nPointsPerPulse} {len(driver.time)}")
                text_spec = []
                freq = []
                if (self.widgets.spectrum_menu.option_var.get() == "Master"):
                    #arrsize = int(2*driver.nPointsPerPulse-1)
                    arrsize = int(2*len(output[1])-1)
                    samplingFrequency   = 1 / driver.Dt
                    freq = np.linspace(-1/2, 1/2,arrsize)*samplingFrequency
                    norm_pulses =np.zeros(arrsize, dtype = complex)
                    for i in range(0,int(driver.nPointsPerPulse)):
                        norm_pulses[i] = output[1][i]
                    
                    Spectrum = abs(np.fft.fftshift( np.fft.fft( norm_pulses)  ) )**2
                    text_spec = np.c_[np.stack([freq, Spectrum], axis = 1)]
                else:
                    if(self.slavepulseid > 0):
                        #arrsize = int(2*driverS.nPointsPerPulse-1)
                        samplingFrequency   = 1 / driver.Dt
                        arrsize = int(2*len(outputs[1])-1)
                        freq = np.linspace(-1/2, 1/2,arrsize)*samplingFrequency
                        norm_pulses =np.zeros(arrsize, dtype = complex)
                        for i in range(0,int(driverS.nPointsPerPulse)):
                            norm_pulses[i] = outputs[1][i]
                        
                        Spectrum = abs(np.fft.fftshift( np.fft.fft( norm_pulses)  ) )**2 
                        text_spec = np.c_[np.stack([freq, Spectrum], axis = 1)]
                
                if (len(freq) > 0): 
                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))
                    axs.set_title("Spectrum")
                    axs.plot(freq,Spectrum)
                    axs.set_xlim([-1000, 1000])
                    axs.set_yscale('log')                
                    axs.grid(which = 'both' , linewidth =1)
                    axs.set(xlabel=r'$Ghz$', ylabel="")
                    fig.tight_layout()
                    plt.show()  
                
                
                np.savetxt('results/last_spectrum_calculation_regular.txt', text_spec)
                
            if (self.widgets.spectrum_on.var.get() == True and self.widgets.spectrumall_on.var.get() == False):
                text_spec =[]
                if (self.widgets.spectrum_menu.option_var.get() == "Master"):
                    Spectrum = wdm.laser_spectrum(driver,output[1],detector_outM,phiout,0,0,driver.nPointsPerPulse,int(pulse1),int(pulse2),0)
                    
                else:
                    Spectrum = wdm.laser_spectrum(driverS,outputs[1],detector_outS,phiSout,0,0,driverS.nPointsPerPulse,int(pulse1),int(pulse2),0)
                          
                if (self.widgets.spectrum2_on.var.get()):
                    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(16, 9))
                    axs[0].set_title("First Pulse")
                    axs[1].set_title("Second Pulse")
                    axs[0].plot(Spectrum[2],Spectrum[0])
                    axs[1].plot(Spectrum[2], Spectrum[1])
                    axs[0].set_xlim([-1000, 1000])
                    axs[1].set_xlim([-1000, 1000])
                    axs[0].set_yscale('log')                
                    axs[1].set_yscale('log')
                    axs[0].grid(which = 'both' , linewidth =1)
                    axs[1].grid(which = 'both' , linewidth =1)
                    axs[1].set(xlabel=r'$Ghz$', ylabel="")
                    fig.tight_layout()
                    plt.show()
                    text_spec = np.c_[np.stack([Spectrum[2], Spectrum[1], Spectrum[0]], axis = 1)]
                else:
                    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))
                    axs.set_title("Spectrum")
                    axs.plot(Spectrum[2],Spectrum[0])
                    axs.set_xlim([-1000, 1000])
                    axs.set_yscale('log')                
                    axs.grid(which = 'both' , linewidth =1)
                    axs.set(xlabel=r'$Ghz$', ylabel="")
                    fig.tight_layout()
                    plt.show() 
                    text_spec = np.c_[np.stack([Spectrum[2], Spectrum[0]], axis = 1)]
                np.savetxt('results/last_spectrum_calculation_regular.txt', text_spec)

           
            print('###################################')

            if (self.widgets.graphic_output.var.get()):
                color = ["tab:orange","tab:green","tab:red","tab:blue"]
                left  = 0.125  # the left side of the subplots of the figure
                right = 0.9    # the right side of the subplots of the figure
                bottom = 0.1   # the bottom of the subplots of the figure
                top =0.95      # the top of the subplots of the figure
                wspace = 0   # the amount of width reserved for blank space between subplots
                hspace = 1   # the amount of height reserved for white space between subplots
            
                     
                fig, axs = plt.subplots(nrows=len(arbplotdata), ncols=1, figsize=(16, 9))
                fig.subplots_adjust(left =left , bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
                
                if(len(arbplotdata) > 1):
                    for i in range(0,len(arbplotdata)): 
                        #if( i == len(arbplotdata)-1):
                        if(arbname[i].find('Interference')!= -1):
                            shift = int(laser.delay / laser.Dt)
                            axs[i].plot(maxtime[0:len(arbplotdata[i])-shift],arbplotdata[i][0:len(arbplotdata[i])-shift],color[i%len(color)])
                        else:
                            axs[i].plot(maxtime,arbplotdata[i],color[i%len(color)])
                        axs[i].set_title(arbname[i])
                        axs[i].yaxis.major.formatter._useMathText = True
                       # print(str(name[i]) + " " + str(len(y[i])) +" " +  str(len(x[0:len(y[i])])))                        
                        if ( i == len(arbplotdata)-1):
                            axs[i].set(xlabel='T ns', ylabel="")
                else:                 
                     axs.plot(maxtime,arbplotdata[0])#,"go",markersize = 3)
                     axs.set_title(arbname)
                     axs.yaxis.major.formatter._useMathText = True
                     axs.grid(which = 'both' , linewidth =1)
                fig.tight_layout()
                plt.show()
                
            
            toolbar = NavigationToolbar2Tk(self.modcanvas, self.widgets.pulseframe_mod, pack_toolbar=False)
            toolbar.update()
            
            self.parent.update()
            inter_delay = 0
            if(dm[6] == 1 or ds[6] == 1):
                inter_delay = int(laser.delay / laser.Dt)
                
            self.widgets.pulsebutton.config(command = Driver.getfig_n(self,toolbar,len(arbplotdata),maxtime,arbplotdata,arbname,arbleg,self.parent, self.modcanvas, self.modfigure,inter_delay,0, 495,540))
            print("Time" , datetime.now() - start)
            
    def stabcalc(self): #perform calculation in the "Stability diagram" frame
        #TODO remove non used parameters 
        # for i in range(0,len(self.widgets.stabinput)):
        #     print(self.widgets.stabinput[i].var.get())    
        ##################################### master driver
        cIb = float(self.data[0])*1e-12
        cIp = float(self.data[1])*1e-12
        crate = float(self.data[3])
        cbandwidth = float(self.data[4])
        cforder = int(self.data[5])
        cfp = float(self.data[6])
        cw = float(self.data[7])        
        clctime = float(self.data[8])
        ##################################### slave driver
        slIb = float(self.slavedata[0])*1e-12
        slIp = float(self.slavedata[1])*1e-12
        slrate =crate                      #slave has the same repetition rate
        slbandwidth = float(self.slavedata[4])
        slorder = int(self.slavedata[5])
        slfp = float(self.slavedata[6])
        slw = float(self.slavedata[7])        
        slctime = float(self.slavedata[8])
        slnpulses = float(self.slavedata[2])
        
        npulses = 0
        arg = []
        slavearg = []

        for i in range(self.datalstart,self.datalend):
            arg.append(float(self.data[i]))
        for i in range(self.datalstart,self.datalend):
            slavearg.append(float(self.slavedata[i]))
            
        arg[25] = float(self.intdata[0])   #interfemoremeter delay
        arg[26] = float(self.intdata[1])   #theta 1
        arg[27] = float(self.intdata[2])   #theta 2
        slavearg[25] = arg[25]
        slavearg[26] = arg[26]
        slavearg[27] = arg[27]
        
        #arg[26] = float(self.intdata[1])   #phase modulation
        #arg[27] = float(self.intdata[2])   #phase modulation 2
        if ( self.data[2] != ""):           #checking pulse number
            npulses = float(self.data[2] )
        #if (npulses_sl > 0):
        #    npulses = npulses_sl
                                           #display time after calculation

        sldriver =  Driver(slIp,slIb,slfp,slw,slctime,slnpulses,slrate) #slave driver
                                                             #initiating laser
        laserS = Laser(cIp,sldriver,slavearg)
        stabvalues = self.stabdata
        output = laserS.stability_diagram(stabvalues[0], stabvalues[1], float(stabvalues[2]), float(stabvalues[3]), float(stabvalues[4]),int(stabvalues[5]))
                
        text_out =[]
        for i in range(0,len(output)):
            text_out.append(output[i])
        self.tnp = np.c_[np.stack(text_out, axis = 1)]  
        np.savetxt('results/last_stab_calculation.txt', self.tnp)
            
        toolbar = NavigationToolbar2Tk(self.stabcanvas, self.widgets.pulseframe_stab, pack_toolbar=False) #a toolbar for the diagram frame
        toolbar.update()
        self.widgets.stabbutton.config(command = Driver.getfig_n(self,toolbar,1,output[0],output[1],"Stability diagram",r"Detuning frequency $, GHz$",self.parent, self.stabcanvas, self.stabfigure,0,self.widgets.pulseframe_stab,600,500))

        if (self.widgets.graphic_output.var.get()): #create a plot in a separate window
              fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))
              left  = 0.125  # the left side of the subplots of the figure
              right = 0.9    # the right side of the subplots of the figure
              bottom = 0.1   # the bottom of the subplots of the figure
              top =0.95      # the top of the subplots of the figure
              wspace = 0   # the amount of width reserved for blank space between subplots
              hspace = 1.2   # the amount of height reserved for white space between subplots
    
              fig.subplots_adjust(left =left , bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
              axs.plot(output[0],output[1],"go")
              axs.set_title("Stability diagram")
              axs.yaxis.major.formatter._useMathText = True
              axs.set(xlabel=r'Injection Ratio, $dB$', ylabel=r"Detuning frequency,$ GHz$")
              fig.tight_layout()
              plt.show()
        #root.update()
        
    def arbpulseadd(self): #create arbitrary pulse of master laser

        if (self.pulseid > self.maxarbpulse):
            self.pulseid = self.maxarbpulse
            return 0
        else:
            if (os.path.exists("config/default_config.json")):
                    f = open("config/default_config.json")
                    data = json.load(f)
                    f.close()
                    self.Arbpulses[self.pulseid].var[0].set(data["arbpulse_settings"]["Ibias"])
                    self.Arbpulses[self.pulseid].var[1].set(data["arbpulse_settings"]["Iptp"])
                    self.Arbpulses[self.pulseid].var[2].set(data["arbpulse_settings"]["pulse_width"])
                    self.Arbpulses[self.pulseid].var[3].set(data["arbpulse_settings"]["gap_current"])
                    self.Arbpulses[self.pulseid].var[4].set(data["arbpulse_settings"]["gap_width"])
                    self.Arbpulses[self.pulseid].var[5].set(data["arbpulse_settings"]["gap_position"])
                    self.Arbpulses[self.pulseid].var[6].set(data["arbpulse_settings"]["gap_repetition_rate"])
                    self.Arbpulses[self.pulseid].var[7].set(data["arbpulse_settings"]["time_shift"])
                    self.Arbpulses[self.pulseid].var[8].set(data["arbpulse_settings"]["integration_rate"])
                    self.Arbpulses[self.pulseid].var[9].set(data["arbpulse_settings"]["incline"])
                    
            for i in range(0,10):
                self.arb_cache[i][self.pulseid] = self.Arbpulses[self.pulseid].var[i].get()

            self.Arbpulses[self.pulseid].put(0+110*(self.pulseid),0,20,self.pulseid,self.widgets.pulseframe_mod_aux)

            
            self.pulseid = self.pulseid + 1
            self.widgets.pulseframe_mod_aux.configure(width = max(self.pulseid,self.slavepulseid)*136)
            self.widgets.mod_frame_canvas.configure(scrollregion = self.widgets.mod_frame_canvas.bbox("all"))
            self.widgets.pulseidlbl.configure( text = f"Number of master pulses : {self.pulseid}")

        '''    
        for i in range(0,70):
            #self.Arbpulses[self.pulseid].put(1140+35*self.pulseid,200,20,self.pulseid)
            self.Arbpulses[self.pulseid].var[2].set("0.8")
            self.Arbpulses[self.pulseid].var[3].set("6")
            self.Arbpulses[self.pulseid].var[4].set("0.2")
            self.Arbpulses[self.pulseid].var[5].set("0.6")
            self.pulseid = self.pulseid + 1
            print(self.pulseid)
        '''    
    def slavearbpulseadd(self): #create arbitrary pulse of slave laser
        if (self.slavepulseid > self.maxarbpulse):
            return 0
        else:
   
            if (os.path.exists("config/default_config.json")):
                f = open("config/default_config.json")
                data = json.load(f)
                f.close()
                self.slaveArbpulses[self.slavepulseid].var[0].set(data["slavearbpulse_settings"]["Ibias"])
                self.slaveArbpulses[self.slavepulseid].var[1].set(data["slavearbpulse_settings"]["Iptp"])
                self.slaveArbpulses[self.slavepulseid].var[2].set(data["slavearbpulse_settings"]["pulse_width"])
                self.slaveArbpulses[self.slavepulseid].var[3].set(data["slavearbpulse_settings"]["gap_current"])
                self.slaveArbpulses[self.slavepulseid].var[4].set(data["slavearbpulse_settings"]["gap_width"])
                self.slaveArbpulses[self.slavepulseid].var[5].set(data["slavearbpulse_settings"]["gap_position"])
                self.slaveArbpulses[self.slavepulseid].var[6].set(data["slavearbpulse_settings"]["gap_repetition_rate"])
                self.slaveArbpulses[self.slavepulseid].var[7].set(data["slavearbpulse_settings"]["time_shift"])
                self.slaveArbpulses[self.slavepulseid].var[8].set(data["slavearbpulse_settings"]["integration_rate"])
                self.slaveArbpulses[self.slavepulseid].var[9].set(data["slavearbpulse_settings"]["incline"])
            for i in range(0,10):
                self.slavearb_cache[i][self.slavepulseid] = self.slaveArbpulses[self.slavepulseid].var[i].get()        
               
            self.slaveArbpulses[self.slavepulseid].put(0+110*self.slavepulseid,480,470,self.slavepulseid,self.widgets.pulseframe_mod_aux)
            #self.slaveArbpulses[self.slavepulseid].var[2].set("0.3")
            #self.slaveArbpulses[self.slavepulseid].var[1].set("20")
            #self.slaveArbpulses[self.slavepulseid].var[6].set("1.25")
            
            self.slaveArbpulses[self.slavepulseid].entry[8].configure(state = "disabled")

            self.slaveArbpulses[self.slavepulseid].entry[8].place_forget()
            self.slaveArbpulses[self.slavepulseid].lbl[8].place_forget()
            
            self.slaveArbpulses[self.slavepulseid].entry[9].place(x = 110*self.slavepulseid, y = 485 +47*8 + 22)
            self.slaveArbpulses[self.slavepulseid].lbl[9].place(x = 110*self.slavepulseid, y = 485 +47*8)
            self.slavepulseid = self.slavepulseid +1
            self.widgets.pulseframe_mod_aux.configure(width = max(self.pulseid,self.slavepulseid)*136)
            self.widgets.mod_frame_canvas.configure(scrollregion = self.widgets.mod_frame_canvas.bbox("all"))
            
           
            self.widgets.slavepulseidlbl.configure(text = f"Number of slave pulses : {self.slavepulseid}")
        
        '''    
        for i in range(0,70):
            #self.slaveArbpulses[self.slavepulseid].put(1140+35*self.slavepulseid,200,470,self.slavepulseid)
            if(self.slavepulseid == 0):
                self.slaveArbpulses[self.slavepulseid].var[7].set("0.2")
            self.slaveArbpulses[self.slavepulseid].var[2].set("0.1")
            self.slaveArbpulses[self.slavepulseid].var[1].set("20")
            self.slavepulseid = self.slavepulseid +1
            print (self.slavepulseid)
        '''    

    def arbpulseclear(self):
        self.pulseid = self.pulseid - 1 
        if (self.pulseid < 0):
            self.pulseid = 0
            return 0
        else:
            #self.pulseid = self.pulseid - 1  
            self.Arbpulses[self.pulseid].hide()
            self.widgets.pulseframe_mod_aux.configure(width = max(self.pulseid,self.slavepulseid)*150)
            self.widgets.pulseidlbl.configure(text = f"Number of master pulses : {self.pulseid}")
    def slavearbpulseclear(self):
        self.slavepulseid = self.slavepulseid - 1 
        if (self.slavepulseid < 0):
            self.slavepulseid = 0
            return 0
        else:
            #self.slavepulseid = self.slavepulseid - 1  
            self.slaveArbpulses[self.slavepulseid].hide()
            self.widgets.pulseframe_mod_aux.configure(width = max(self.pulseid,self.slavepulseid)*150)
            self.widgets.slavepulseidlbl.configure(text = f"Number of slave pulses : {self.slavepulseid}")

    def callback(self,e): #test function, not used
        self.mx_p = self.mx
        self.my_p = self.my
        
        self.mx = e.x
        self.my = e.y
        if(abs(self.mx - self.mx_p) > 20 and abs(self.my - self.my_p) > 20  ):
            print("value reset")
        # x= e.x
        # y= e.y
        print("Pointer is currently at %d, %d" %(self.mx,self.my))
    def buttonPressed( self, event ): #test function, not used
        print( "Pressed at [ " + str( event.x ) + 
                              ", " + str( event.y ) + " ]" )

    def savelaser(self, event = None): #main function that saves values in all fields
        for i in range(0,self.intstart):
            if (self.entry[i].get() != ''):
                self.data[i] = self.entry[i].get()
                
        
            if (self.slaveentry[i].get() != ''):
                self.slavedata[i] = self.slaveentry[i].get()
        
        for i in range(self.intstart,self.datalend):
            self.data[i] = self.intdata[i- self.intstart]
            
        for i in range(0,3):
            if (self.dentry[i].get() != ''):
                self.wdmdata[i] = self.dentry[i].get() 
        for i in range(0,3):
            if (self.slavedentry[i].get() != ''):
                self.slavewdmdata[i] = self.slavedentry[i].get() 
        for i in range(0,3):
            if (self.intentry[i].get() != ''):
                self.intdata[i] = self.intentry[i].get()
        for i in range(0,3):
            if (self.photentry[i].get() != ''):
                self.photdata[i] = self.photentry[i].get()
        
        for i in range(0,self.datadend):
                self.data[i] = self.widgets.master_entry[i].spinboxm.get()
                self.slavedata[i] = self.widgets.slave_entry[i].spinboxm.get()
                
        for i in range(0,len(self.widgets.stabinput)):
            self.stabdata[i] = self.widgets.stabinput[i].var.get()     
        
        for i in range(0, len(self.widgets.transfer_range)):
            self.transdata[i] = self.widgets.transfer_range[i].var.get()
        for j in range(0,self.pulseid):
            for i in range(0,10):
                    self.arb_cache[i][j] = self.Arbpulses[j].var[i].get()
                    
        for j in range(0,self.slavepulseid):
            for i in range(0,10):            
                    self.slavearb_cache[i][j] = self.slaveArbpulses[j].var[i].get()
        # for i in range(0,self.datalend):
        #     print(self.data[i])
        print("parameters saved")        

          
    def driver_reset(self,event): #function that resets value in driver fields

        for i in range(0,self.datadend):
            self.widgets.master_entry[i].spinboxm.set(self.data[i])
            self.widgets.slave_entry[i].spinboxm.set(self.slavedata[i])
            self.widgets.master_entry[i].spinboxm.update()
            self.widgets.slave_entry[i].spinboxm.update()    
        print("driver reset ")
        
    def laser_reset(self,event):                       
        for i in range(9, 34):
            self.entry[i].set(self.data[i])
            self.entry[i].update()
        for i in range(9, 34):
            self.slaveentry[i].set(self.slavedata[i])
            self.slaveentry[i].update()
        print("laser reset ")
    def arbpulse_reset(self,event):
        for i in range(0,self.pulseid): 
            for j in range(0,10):
                self.Arbpulses[i].entry[j].set( self.arb_cache[j][i])
                #print( self.arb_cache[j][i])
                self.Arbpulses[i].entry[j].update()
        print("arbpulse reset ")
    def slavearbpulse_reset(self,event):
        for i in range(0,self.slavepulseid):   
            for j in range(0,10): 
                self.slaveArbpulses[i].entry[j].set(self.slavearb_cache[j][i])   
                
                self.slaveArbpulses[i].entry[j].update()  
        print("slavearbpulse reset ")
    def input_reset(self,event):
        
        print("input reset")        
        
        self.intvar[0].set(self.intdata[0])
        self.intvar[1].set(self.intdata[1])
        self.intvar[2].set(self.intdata[2])
        
        for i in range(0,3):
            self.dentry[i].set(self.wdmdata[i])
            self.dentry[i].update()
        for i in range(0,3):
            self.slavedentry[i].set(self.slavewdmdata[i])
            self.slavedentry[i].update()
        for i in range(0,3):
            
            self.intentry[i].update()
        for i in range(0,3):
            self.photentry[i].set(self.photdata[i])
            self.photentry[i].update()
        for i in range(0,self.pulseid+1):
            for j in range(0,10):
                self.Arbpulses[i].entry[j].update()
        for i in range(0,self.slavepulseid+1):
            for j in range(0,10):
                self.slaveArbpulses[i].entry[j].update()
                
        for i in range(0,len( self.widgets.stabinput)):
            self.widgets.stabinput[i].spinboxm.set(self.stabdata[i])
            self.widgets.stabinput[i].spinboxm.update()
        for i in range(len(self.widgets.transfer_range)):
            self.widgets.transfer_range[i].var.set(self.transdata[i])
            self.widgets.transfer_range[i].spinboxm.update()
        
            #self.parent.response_wp_1.update()
            #self.parent.response_wp_1.config( textvariable = '0')
            #f.update()
            #print("focus out" )            
    def config_train(self):
        temp_arr = conf_pulsetrain(self.pulseid,self.slavepulseid,self.arb_cache,self.slavearb_cache,"trainconfiguration.json")
        filename = filedialog.asksaveasfilename()

        if filename != None:
            if filename != '':
                with open(str(filename),"w") as write_file:
                    json.dump(temp_arr,write_file,indent=4)
        
        #TODO the train list only updates on launch
        # dir_path = r'config/trains'
        # res = []
        # res.append("New train")
        # for path in os.listdir(dir_path):
        # # check if current path is a file
        #     if os.path.isfile(os.path.join(dir_path, path)):
        #         res.append(path)
        # current_var = self.widgets.train_menu.option_var.get()
        # self.widgets.train_menu = Optionbox(self.widgets.pulseframe_mod,res , change = self.train_change, tip = "Select saved train from \ config \ trains folder\nTo clean or make a new train, select - New train ")
        # self.widgets.train_menu.option_var.set(current_var)
        # self.widgets.train_menu.calc_menu.place (x=0 + 260, y=620)
        
    def loadconfig_train(self, file = None):
        if (file == None):
            filename = filedialog.askopenfilename()
        else:
            filename = file
        data =[]
        
        if filename != '':
            f = open(filename)
            data = json.load(f)
            f.close()
            
            pulse_num = data[0]["pulse_num"]
            slavepulse_num = data[0]["slavepulse_num"]
            self.widgets.train_menu.option_var.set(Path(filename).name)
            dir_path = r'config/trains'
            res = [] #updated list of trains
            res.append("New train")
            for path in os.listdir(dir_path):
            # check if current path is a file
                if os.path.isfile(os.path.join(dir_path, path)):
                    res.append(path) 
            self.widgets.train_menu.options = res
            #self.widgets.train_menu.calc_menu.update()
            #self.initUI 
            #TODO the train list only updates on launch
            if (self.pulseid-pulse_num > 0):      #if there are more pulses right now compared to train clean the rest
                for i in range(0, self.pulseid-pulse_num):
                    #print("delete pulse")
                    self.arbpulseclear()
            if (self.slavepulseid-slavepulse_num > 0): # else add more pulses
                for i in range(0, self.slavepulseid-slavepulse_num):
                    print("delete slavepulse")
                    #self.slavearbpulseclear()
            if (self.pulseid-pulse_num < 0):
                for i in range(0,pulse_num - self.pulseid):
                    self.arbpulseadd()        
            
            for i in range(0,pulse_num):         
                self.Arbpulses[i].var[0].set(data[1][i]["arbpulse_settings"]["Ibias"])
                self.Arbpulses[i].var[1].set(data[1][i]["arbpulse_settings"]["Iptp"])
                self.Arbpulses[i].var[2].set(data[1][i]["arbpulse_settings"]["pulse_width"])
                self.Arbpulses[i].var[3].set(data[1][i]["arbpulse_settings"]["gap_current"])
                self.Arbpulses[i].var[4].set(data[1][i]["arbpulse_settings"]["gap_width"])
                self.Arbpulses[i].var[5].set(data[1][i]["arbpulse_settings"]["gap_position"])
                self.Arbpulses[i].var[6].set(data[1][i]["arbpulse_settings"]["gap_repetition_rate"])
                self.Arbpulses[i].var[7].set(data[1][i]["arbpulse_settings"]["time_shift"])
                self.Arbpulses[i].var[8].set(data[1][i]["arbpulse_settings"]["integration_rate"])
                self.Arbpulses[i].var[9].set(data[1][i]["arbpulse_settings"]["incline"])
        
            if (self.slavepulseid - slavepulse_num < 0):
                for i in range(0, slavepulse_num - self.slavepulseid):
                    self.slavearbpulseadd()

            for i in range(0,slavepulse_num):
                self.slaveArbpulses[i].var[0].set(data[2][i]["slavearbpulse_settings"]["Ibias"])
                self.slaveArbpulses[i].var[1].set(data[2][i]["slavearbpulse_settings"]["Iptp"])
                self.slaveArbpulses[i].var[2].set(data[2][i]["slavearbpulse_settings"]["pulse_width"])
                self.slaveArbpulses[i].var[3].set(data[2][i]["slavearbpulse_settings"]["gap_current"])
                self.slaveArbpulses[i].var[4].set(data[2][i]["slavearbpulse_settings"]["gap_width"])
                self.slaveArbpulses[i].var[5].set(data[2][i]["slavearbpulse_settings"]["gap_position"])
                self.slaveArbpulses[i].var[6].set(data[2][i]["slavearbpulse_settings"]["gap_repetition_rate"])
                self.slaveArbpulses[i].var[7].set(data[2][i]["slavearbpulse_settings"]["time_shift"])
                self.slaveArbpulses[i].var[8].set(data[2][i]["slavearbpulse_settings"]["integration_rate"])
                self.slaveArbpulses[i].var[9].set(data[2][i]["slavearbpulse_settings"]["incline"])
            
            for i in range(0,pulse_num):
                for j in range(0,10):
                   self.slaveArbpulses[i].entry[j].update() 
                   
            for i in range(0,slavepulse_num):
                for j in range(0,10):
                   self.slaveArbpulses[i].entry[j].update() 

            
    def config_b(self):
        
            for  i in range(0,self.intstart):
                if (self.entry[i].get() != ''):
                    self.data[i] = self.entry[i].get()
                if (self.slaveentry[i].get() != ''):
                    self.slavedata[i] = self.slaveentry[i].get()
            for i in range(self.intstart,self.datalend):
                self.data[i] = self.intdata[i- self.intstart]
                
            for i in range(0,3):
                if (self.dentry[i].get() != ''):
                    self.wdmdata[i] = self.dentry[i].get() 
            for i in range(0,3):
                if (self.slavedentry[i].get() != ''):
                    self.slavewdmdata[i] = self.slavedentry[i].get() 
            for i in range(0,3):
                if (self.intentry[i].get() != ''):
                    self.intdata[i] = self.intentry[i].get()
            for i in range(0,3):
                if (self.photentry[i].get() != ''):
                    self.photdata[i] = self.photentry[i].get()
            # for i in range(0,7):
            #     if (self.arbentry[i].get() != ''):
            #         self.arbdata[i] = self.arbentry[i].get() 
                    
            for i in range(0,len(self.widgets.stabinput)):
                self.stabdata[i] = self.widgets.stabinput[i].var.get() 
            
            arbdata = np.zeros(10)
            slavearbdata = np.zeros(10)
            for i in range(0,10): 
                arbdata[i] = self.Arbpulses[0].entry[i].get()
            for i in range(0,10): 
                slavearbdata[i] = self.slaveArbpulses[0].entry[i].get()
            for i in range(0,len(self.widgets.transfer_range)):
                self.transdata[i] = self.widgets.transfer_range[i].var.get()
            #trans_arr = [self.widgets.transfer_range[0].var.get(),self.widgets.transfer_range[1].var.get()]#range for transfer func
            opt_arr = [self.widgets.f_var.get(),  self.widgets.wdmvar1.get(), self.widgets.wdmvar2.get(),self.widgets.calc_var.get()]
            temp_arr = config(self.data,self.wdmdata,self.slavedata,self.slavewdmdata,self.photdata,arbdata,slavearbdata,self.stabdata,self.transdata,opt_arr,"configuration.json")
            filename = filedialog.asksaveasfilename()

            if filename != None:
                if filename != '':
                    with open(str(filename),"w") as write_file:
                        json.dump(temp_arr,write_file)
            else:
                tk_messagebox.showinfo("Write filename")
    def config_setup(self):         #getting values from the configurations
            filename = filedialog.askopenfilename()
            if filename != '':
                f = open(filename)
                data = json.load(f)
                with open("config/default_config.json","w") as write_file:
                    json.dump(data,write_file)
                f.close()
            
                if (os.path.exists("config/default_config.json")):
                    f = open("config/default_config.json")
                    data = json.load(f)
                    self.myvar[0].set(data["driver_settings"]["Ibias"]*1e12)
                    self.myvar[1].set(data["driver_settings"]["Iptp"]*1e12)
                    self.myvar[2].set(
                        data["driver_settings"]["number of Pulses"])
                    self.myvar[3].set(
                        data["driver_settings"]["integration rate"])
                    self.myvar[4].set(
                        data["driver_settings"]["filter bandwidth"])
                    self.myvar[5].set(data["driver_settings"]["filter order"])
                    self.myvar[6].set(
                        data["driver_settings"]["pulse repetiotion rate"])
                    self.myvar[7].set(data["driver_settings"]["pulse width"])
                    self.myvar[8].set(data["driver_settings"]["decay time"])
                    
                    self.slavevar[0].set(
                        data["slave_driver_settings"]["Ibias"]*1e12)
                    self.slavevar[1].set(
                        data["slave_driver_settings"]["Iptp"]*1e12)
                    self.slavevar[2].set(
                        data["slave_driver_settings"]["number of Pulses"])
                    self.slavevar[3].set(
                        data["slave_driver_settings"]["integration rate"])
                    self.slavevar[4].set(
                        data["slave_driver_settings"]["filter bandwidth"])
                    self.slavevar[5].set(
                        data["slave_driver_settings"]["filter order"])
                    self.slavevar[6].set(
                        data["slave_driver_settings"]["pulse repetiotion rate"])
                    self.slavevar[7].set(
                        data["slave_driver_settings"]["pulse width"])
                    self.slavevar[8].set(
                        data["slave_driver_settings"]["decay time"])
                    
                    self.myvar[9].set(data["laser_settings"]["eta"])
                    self.myvar[10].set(data["laser_settings"]["t_ph"])
                    self.myvar[11].set(data["laser_settings"]["N_transparency"])
                    self.myvar[12].set(data["laser_settings"]["N_threshhold"])
                    self.myvar[13].set(data["laser_settings"]["chi"])
                    self.myvar[14].set(data["laser_settings"]["C_sp"])
                    self.myvar[15].set(data["laser_settings"]["Gamma"])
                    self.myvar[16].set(data["laser_settings"]["tau_e"])
                    self.myvar[17].set(data["laser_settings"]["alpha"])
                    self.myvar[18].set(data["laser_settings"]["rh"])
                    self.myvar[19].set(data["laser_settings"]["th"])
                    self.myvar[20].set(data["laser_settings"]["lambda"])
                    self.myvar[21].set(data["laser_settings"]["kap_om"])
                    self.myvar[22].set(data["laser_settings"]["t1"])
                    self.myvar[23].set(data["laser_settings"]["w_M"])
                    self.myvar[24].set(data["laser_settings"]["om_th"])
                    self.myvar[25].set(data["laser_settings"]["Ll"])
                    self.myvar[26].set(data["laser_settings"]["kap_inj"])
                    self.myvar[27].set(data["laser_settings"]["r1"])
                    self.myvar[28].set(data["laser_settings"]["r2"])
                    self.myvar[29].set(data["laser_settings"]["rhs"])
                    self.myvar[30].set(data["laser_settings"]["k12"])
                    self.myvar[31].set(data["laser_settings"]["tauj"])
                    self.myvar[32].set(data["laser_settings"]["lambdas"])
                    self.myvar[33].set(data["laser_settings"]["Cj"])
                    self.myvar[34].set(data["laser_settings"]["int_delay"])
                    self.myvar[35].set(data["laser_settings"]["phi1"])
                    self.myvar[36].set(data["laser_settings"]["phi2"])
                    
                    self.slavevar[9].set(data["slave_laser_settings"]["eta"])
                    self.slavevar[10].set(data["slave_laser_settings"]["t_ph"])
                    self.slavevar[11].set(
                        data["slave_laser_settings"]["N_transparency"])
                    self.slavevar[12].set(
                        data["slave_laser_settings"]["N_threshhold"])
                    self.slavevar[13].set(data["slave_laser_settings"]["chi"])
                    self.slavevar[14].set(data["slave_laser_settings"]["C_sp"])
                    self.slavevar[15].set(data["slave_laser_settings"]["Gamma"])
                    self.slavevar[16].set(data["slave_laser_settings"]["tau_e"])
                    self.slavevar[17].set(data["slave_laser_settings"]["alpha"])
                    self.slavevar[18].set(data["slave_laser_settings"]["rh"])
                    self.slavevar[19].set(data["slave_laser_settings"]["th"])
                    self.slavevar[20].set(
                        data["slave_laser_settings"]["lambda"])
                    self.slavevar[21].set(
                        data["slave_laser_settings"]["kap_om"])
                    self.slavevar[22].set(data["slave_laser_settings"]["t1"])
                    self.slavevar[23].set(data["slave_laser_settings"]["w_M"])
                    self.slavevar[24].set(data["slave_laser_settings"]["om_th"])
                    self.slavevar[25].set(data["slave_laser_settings"]["Ll"])
                    self.slavevar[26].set(
                        data["slave_laser_settings"]["kap_inj"])
                    self.slavevar[27].set(data["slave_laser_settings"]["r1"])
                    self.slavevar[28].set(data["slave_laser_settings"]["r2"])
                    self.slavevar[29].set(data["slave_laser_settings"]["rhs"])
                    self.slavevar[30].set(data["slave_laser_settings"]["k12"])
                    self.slavevar[31].set(data["slave_laser_settings"]["tauj"])
                    self.slavevar[32].set(
                        data["slave_laser_settings"]["lambdas"])
                    self.slavevar[33].set(data["slave_laser_settings"]["Cj"])
                    self.slavevar[34].set(
                        data["slave_laser_settings"]["int_delay"])
                    self.myvar[35].set(data["laser_settings"]["phi1"])
                    self.myvar[36].set(data["laser_settings"]["phi2"])
                    
                    self.wdmvar[0].set(data["wdm_settings"]["filter bandwidth"])
                    self.wdmvar[1].set(data["wdm_settings"]["filter order"])
                    self.wdmvar[2].set(data["wdm_settings"]["filter shift"])
                    
                    self.slavewdmwar[0].set(
                        data["slave_wdm_settings"]["filter bandwidth"])
                    self.slavewdmwar[1].set(
                        data["slave_wdm_settings"]["filter order"])
                    self.slavewdmwar[2].set(
                        data["slave_wdm_settings"]["filter shift"]) 
                    self.intvar[0].set(self.myvar[34].get())
                    self.intvar[1].set(self.myvar[35].get())
                    self.intvar[2].set(self.myvar[36].get())
                    
                    self.widgets.stabinput[0].var.set(data["diagram_settings"]["Qinj_min"])
                    self.widgets.stabinput[1].var.set(data["diagram_settings"]["Qinj_max"])
                    self.widgets.stabinput[2].var.set(data["diagram_settings"]["phi_min"])
                    self.widgets.stabinput[3].var.set(data["diagram_settings"]["phi_max"])
                    self.widgets.stabinput[4].var.set(data["diagram_settings"]["sol_point"])
                    self.widgets.stabinput[5].var.set(data["diagram_settings"]["pointsnum"])
                   
                    for i in range(0,self.pulseid+1):
                        self.Arbpulses[i].var[0].set(data["arbpulse_settings"]["Ibias"])
                        self.Arbpulses[i].var[1].set(data["arbpulse_settings"]["Iptp"])
                        self.Arbpulses[i].var[2].set(data["arbpulse_settings"]["pulse_width"])
                        self.Arbpulses[i].var[3].set(data["arbpulse_settings"]["gap_current"])
                        self.Arbpulses[i].var[4].set(data["arbpulse_settings"]["gap_width"])
                        self.Arbpulses[i].var[5].set(data["arbpulse_settings"]["gap_position"])
                        self.Arbpulses[i].var[6].set(data["arbpulse_settings"]["gap_repetition_rate"])
                        self.Arbpulses[i].var[7].set(data["arbpulse_settings"]["time_shift"])
                        self.Arbpulses[i].var[8].set(data["arbpulse_settings"]["integration_rate"])
                        self.Arbpulses[i].var[9].set(data["arbpulse_settings"]["incline"])
                    for i in range(0,self.slavepulseid+1):
                        
                        self.slaveArbpulses[i].var[0].set(data["slavearbpulse_settings"]["Ibias"])
                        self.slaveArbpulses[i].var[1].set(data["slavearbpulse_settings"]["Iptp"])
                        self.slaveArbpulses[i].var[2].set(data["slavearbpulse_settings"]["pulse_width"])
                        self.slaveArbpulses[i].var[3].set(data["slavearbpulse_settings"]["gap_current"])
                        self.slaveArbpulses[i].var[4].set(data["slavearbpulse_settings"]["gap_width"])
                        self.slaveArbpulses[i].var[5].set(data["slavearbpulse_settings"]["gap_position"])
                        self.slaveArbpulses[i].var[6].set(data["slavearbpulse_settings"]["gap_repetition_rate"])
                        self.slaveArbpulses[i].var[7].set(data["slavearbpulse_settings"]["time_shift"])
                        self.slaveArbpulses[i].var[8].set(data["slavearbpulse_settings"]["integration_rate"])
                        self.slaveArbpulses[i].var[9].set(data["slavearbpulse_settings"]["incline"])
                    
                    self.widgets.transfer_range[0].var.set(data["transfer_settings"]["freq_range_start"])
                    self.widgets.transfer_range[1].var.set(data["transfer_settings"]["freq_range_end"])
                    
                    f.close
                    
                    
                else:
                    tk_messagebox.showinfo(
                        "Error", "configuration file has not been found!")
                    for i in range(0, self.datadend):
                        self.myvar[i].set("0")
                        
                for i in range(0,len( self.widgets.stabinput) ):
                    for i in range(0,len( self.widgets.stabinput)):         
                        self.stabdata[i] = self.widgets.stabinput[i].var.get()
                        self.widgets.stabinput[i].spinboxm.update()
                        
                for i in range(0, self.datadend):
                    self.entry[i].update()
                for i in range(0, self.datadend):
                    if (self.slavelbl[i].cget("text") !="Integration rate"):
                        self.slaveentry[i].update()
                        
                for i in range(0,self.datadend):
                    #print(self.data[i])
                    #self.myvar[i].set(self.data[i])
                    self.widgets.master_entry[i].spinboxm.set(self.myvar[i].get())
                    self.widgets.slave_entry[i].spinboxm.set(self.slavevar[i].get())
                    self.widgets.master_entry[i].spinboxm.update()
                    self.widgets.master_entry[i].spinboxm.update()
                for i in range(9, 34):
                     self.entry[i].update()
                for i in range(9, 34):
                    self.slaveentry[i].update()
                    
                for i in range(0,3):
                    self.dentry[i].update()
                for i in range(0,3):
                    self.slavedentry[i].update()
                for i in range(0,3):
                    self.intentry[i].update()   
                for i in range(0,self.pulseid+1):
                    for j in range(0,10):
                        self.Arbpulses[i].entry[j].update()
                for i in range(0,self.slavepulseid+1):
                    for j in range(0,10):
                        self.slaveArbpulses[i].entry[j].update()
                for i in range(len(self.widgets.transfer_range)):
                    self.transdata[i] = self.widgets.transfer_range[i].var.get()
                    self.widgets.transfer_range[i].spinboxm.update() 
                #print(data["driver_settings"]["Ibias"])
                #f.close
    def train_change(self, *args): #clean pulses field or load a train
        
        if (self.widgets.train_menu.option_var.get() != "New train"):
            filename = r'config/trains/' + str(self.widgets.train_menu.option_var.get())
            self.widgets.train_menu.calc_menu.update()
            self.loadconfig_train(filename)
        if (self.widgets.train_menu.option_var.get() == "New train"):
            for i in range(0,self.pulseid):
                self.arbpulseclear()
            for i in range(0,self.slavepulseid):
                self.slavearbpulseclear()
    def option_changed(self, *args): #displays and hides elements depending on the chosen scheme
        self.widgets.label['text'] = f'Current composition: {self.widgets.option_var.get()}'
        
        if(self.widgets.option_var.get().find('Slave') != -1 ):
            self.slavestate = 1            
            c = 10
            for i in range(0, self.datadend):
                if ( self.widgets.slave_entry[i].label.cget("text") !="Integration rate"):
                    self.widgets.slave_entry[i].label.place(x = 5, y= c + 3)           
                    self.widgets.slave_entry[i].spinboxm.place(x = 120, y= c, width=100)
                    c = c + 35                        
            c = 10
            

            for i in range(9, 18):
                self.slaveentry[i].place(x=90, y=c, width=110)
                self.slavelbl[i].place(x=0, y = c + 3)
                c = c + 35
            c = 10    
            for i in range(18, 26):
                if( (i in self.widgets.unnecessary) != 1 ):
                    self.slaveentry[i].place(x=150, y=c, width=110)
                    self.slavelbl[i].place(x=0, y = c + 3)
                    c = c + 35
            c = 10    
            for i in range(26, 34):
                if( (i in self.widgets.unnecessary) != 1 ):
                    self.slaveentry[i].place(x=370, y=c, width=110)
                    self.slavelbl[i].place(x=320, y = c + 3)
                    c = c + 35
            c = 290            
            for i in range(0, self.widgets.slave_clcoptions):
                self.widgets.s_checkbox_clc[i].checkboxm.place(x = 10 , y = 10 + c)
                c = c + 30

            c = 270
            for i in range(0, self.widgets.slave_pltoptions):
                self.widgets.s_checkbox_plt[i].checkboxm.place(x = 10, y = 10 + c)
                c = c + 30            
            self.widgets.lnotebooks.place(x = 10, y = 10)
        else:                   #Slave laser should be hidden
            self.slavestate = 0
            for i in range(0, self.widgets.slave_clcoptions):
                self.widgets.s_checkbox_clc[i].checkboxm.place_forget()
            for i in range(0, self.widgets.slave_pltoptions):
                self.widgets.s_checkbox_plt[i].checkboxm.place_forget()
            self.widgets.lnotebooks.place_forget()
            for i in range(0, self.datadend):
                if (self.widgets.slave_entry[i].label.cget("text") !="Integration rate"):
                    self.widgets.slave_entry[i].label.place_forget()         
                    self.widgets.slave_entry[i].spinboxm.place_forget()
            for i in range(self.datadend, 34):
                self.slaveentry[i].place_forget()
                self.slavelbl[i].place_forget()
        if(self.widgets.option_var.get().find('Interferometer') != -1 ):    #display interferometer input
            c = self.wdmy   - 20               
            for i in range(0,3):
                self.intentry[i].place(x = self.wdmx+275, y = c + 26, width = 100)
                self.intlbl[i].place(x = self.wdmx+275, y = c + 3  )
                c  = c + 52
        else:
            for i in range(0,3):
                self.intentry[i].place_forget()
                self.intlbl[i].place_forget()
            
        if(self.widgets.option_var.get().find('Master-WDM') != -1 ):  #display master WDM input
            c = self.wdmy - 20
            for i in range(0,3):
                self.dentry[i].place(x = self.wdmx-30, y = c + 27, width = 100)
                self.setuplbl[i].place(x = self.wdmx-30, y = c + 5 )
                c  = c + 50
        else:
            for i in range(0,3):
                self.dentry[i].place_forget()
                self.setuplbl[i].place_forget()
        if(self.widgets.option_var.get().find('Slave-WDM') != -1 ):    #display slave WDM input
            c = self.wdmy - 20
            for i in range(0,3):
                self.slavedentry[i].place(x = self.wdmx+105, y = c + 27, width = 100)
                self.slavewdmlbl[i].place(x = self.wdmx+105, y = c + 5)
                c  = c + 50
        else:
            for i in range(0,3):
                self.slavedentry[i].place_forget()
                self.slavewdmlbl[i].place_forget()       
        imgx = 700
        imgy= 5
        for i in range(0,10):                           #Changing scheme pictures
            self.Imlbl[i].place_forget()
        if (self.widgets.option_var.get() == "Driver-Master"):
            self.Imlbl[0].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-WDM"):
            self.Imlbl[1].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-Interferometer"):
            self.Imlbl[2].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-WDM-Interferometer"):
            self.Imlbl[3].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-Slave"):
            self.Imlbl[4].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-Slave-Interferometer"):
            self.Imlbl[5].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-Slave-WDM"):
            self.Imlbl[6].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-Slave-WDM-Interferometer"):
            self.Imlbl[7].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == "Driver-Master-WDM-Slave-WDM"):
            self.Imlbl[8].place(x = imgx , y =imgy)
        if (self.widgets.option_var.get() == 'Driver-Master-WDM-Slave-WDM-Interferometer'):
            self.Imlbl[9].place(x = imgx , y =imgy)

    def filter_changed(self, *args):
       
        if(self.widgets.f_var.get() == "Bessel" ):
            self.filter_driver  = "Bessel"
        elif(self.widgets.f_var.get() =="Butterworth"):
            self.filter_driver  = "Butterworth"
        if (self.widgets.calc_var.get() == "Runge_Kutta"):
            
            self.widgets.calclabel.config(text = "Calculation scheme (non stochastic)")
           
            #self.widgets.Euler_mod.calc_menu.place_forget()
            #self.widgets.scipym.checkboxm.place(x=760, y=780 + 50)
        elif (self.widgets.calc_var.get() == "Euler"):
            #if (self.widgets.m_checkbox_clc[5].var.get() == 1):
            
            self.widgets.calclabel.config(text = "Calculation scheme (stochastic)")
            #self.widgets.Euler_mod.calc_menu.place( x = 880, y = 780)
            #self.widgets.scipym.checkboxm.place_forget()
            #x=760, y=780
    def clearme(self):       
        self.canvas.get_tk_widget().place_forget() 
    def download_spectrum(self): 
        if (self.widgets.s1.get()):
           filename = filedialog.asksaveasfile(mode='w',defaultextension=".csv")
        else:
            filename = filedialog.asksaveasfile(mode='w',defaultextension=".txt")
        
        if filename is not None:
            if filename != '' :
                if(self.widgets.s1.get()):
                    #np.savetxt(filename , [self.tnp[:,1]],delimiter=',')
                    #np.transpose(self.tnp).tofile(filename, sep =',')
                    DF = pd.DataFrame(self.tmp_spec) 
                    DF.to_csv(filename , sep = ',', index = 0)
                else:
                    np.savetxt(filename , self.tmp_spec )
            else:
                #TODO
                today = date.today()                
                np.savetxt("SIMLAD results: " + str(today),self.tmp_spec)
                
                
    def download(self): 
        if (self.widgets.s1.get()):
           filename = filedialog.asksaveasfile(mode='w',defaultextension=".csv")
        else:
            filename = filedialog.asksaveasfile(mode='w',defaultextension=".txt")
        
        if filename is not None:
            if filename != '' :
                if(self.widgets.s1.get()):
                    #np.savetxt(filename , [self.tnp[:,1]],delimiter=',')
                    #np.transpose(self.tnp).tofile(filename, sep =',')
                    DF = pd.DataFrame(self.tnp) 
                    DF.to_csv(filename , sep = ',', index = 0)
                else:
                    np.savetxt(filename , self.tnp )
            else:
                #TODO
                today = date.today()                
                np.savetxt("SIMLAD results: " + str(today),self.tnp)
                
                
    def initUI(self):
        
        self.parent.title("SIMLAD")

        self.place(x = 0 , y = 0)
        self.widgets = Widget(self)
        #self.widgets.grid(row = 3 , column = 3)
        #self.widgets.pack(side="top", anchor="center", fill="both", expand=True)

    def change(self):
        for i in range(self.datadend,self.datalend):
            if (self.entry[i].winfo_exists()):
                if (self.entry[i].get() != ''):
                    self.data[i] = self.entry[i].get()
        for i in range(0, self.datadend):
            if (self.widgets.master_entry[i].var.get() != ''):
                self.data[i] = self.widgets.master_entry[i].var.get()
                
    def slchange(self):
        for i in range(self.datadend,self.datalend):
            if (self.slaveentry[i].winfo_exists()):
                if (self.slaveentry[i].get() != ''):
                    self.slavedata[i] = self.slaveentry[i].get()
        for i in range(0, self.datadend):
           if (self.widgets.slave_entry[i].var.get() != ''):
               self.slavedata[i] = self.widgets.slave_entry[i].var.get()
                          
    def clear(self):
        for item in  self.canvas.get_tk_widget().find_all():
            self.canvas.get_tk_widget().delete(item)
        self.parent.update()
              


root = themed_tk.ThemedTk() #tk.Tk()#
#root.get_themes()                 # Returns a list of all themes that can be set
#root.set_theme("equilux") 

style = ThemedStyle(root)
style.set_theme("arc") #arc, equilux, radiance

root.geometry("1900x1000+0+0")
app = App(root)
app.pack(side="top", fill="both", expand=True)
root.bind("<Return>", lambda event: app.savelaser(event = event))
root.mainloop()