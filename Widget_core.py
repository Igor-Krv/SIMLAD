import tkinter as tk
import numpy as np
from tkinter import ttk
from Interface_elements_core import Input
from Interface_elements_core import Checkbox
from Interface_elements_core import Optionbox
from tkinter import messagebox as tk_messagebox
from idlelib.tooltip import Hovertip
#from IPython.display import display, Latex
from datetime import datetime
import json
import os.path
from PIL import ImageTk, Image
import os

class Widget(ttk.Frame):
    def __init__(self, parent):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        self.lstate = 0    # laser window state
        self.dstate = 0    # driver window state
        self.cstate = 0    # calculation window state
        self.pstate = 0    # plot window state
        self.setup = 0     # system setup  state
        self.wdmstate = 0  # WDM filter window state
        self.intstate = 0  # inteferometer window state
        
        self.initUI()
        self.initDriver()
        self.initSlaveDriver()
        self.initLaser()
        self.initWDM()
        self.initInt()
        self.initSlaveLaser()
        self.initSlaveWDM()
        self.InitPhot()
        self.initModFrame()
        self.initStabFrame()
    def empty_scroll_command(self, event):
         return "break"   

         
    def initUI(self):
        
        def commands(event,f):
            self.parent.response_wp_1.update()
            self.parent.response_wp_1.config( textvariable = '0')
            f.update()
            print("focus out" )
         

        self.checkbx_plt = []
        self.entry = []
        for i in range(1, 44):
            self.checkbx_plt.append(ttk.Checkbutton())
            self.entry.append(ttk.Spinbox())

        pulseframe_x = 1010
        pulseframe_y = 10
        self.pulseframe = ttk.Notebook(self.parent)
        self.pulseframe.place(x=pulseframe_x, y=pulseframe_y)
        self.pulseframe_reg = tk.Frame(self.pulseframe, width=900, height=960)

        self.pulseframe.add(self.pulseframe_reg, text='Regular Simulations')

        self.button = ttk.Button(self.pulseframe_reg, text='Build plot',
                                 command=self.parent.maincalculation, cursor="hand2")
        self.button.place(x=0, y = 750, width = 150, height = 50)
        
        self.spectrumreg_on = Checkbox(self.pulseframe_reg,"Calculate spectrum of the train")    
        self.spectrumreg_on.checkboxm.place(x= 0, y= 750 + 64)
        
        self.download_spec = ttk.Button(self.pulseframe_reg, text='Download spectrum data',
                               command=self.parent.download_spectrum, cursor="hand2")
        #self.download_spec.place(x = 250, y= 750, width = 250, height = 50)
        
        self.spectrum_menu_reg = Optionbox(self.pulseframe_reg,('Master', 'Slave') )
        self.spectrum_menu_reg.calc_menu.place (x = 0, y = 750 + 128)
        
        self.specpulse_s = Input(self.pulseframe_reg, 0, 30, 1, 1)
        self.specpulse_s.spinboxm.config(width = 10)
        self.specpulse_s.spinboxm.place(x = 100 , y = 750 + 108)
        
        self.specpulse_end = Input(self.pulseframe_reg, 0, 30, 1, 1)
        self.specpulse_end.spinboxm.config(width = 10)
        self.specpulse_end.spinboxm.place(x = 100, y = 750 + 148)
        
        # self.warninglabel = ttk.Label(self.pulseframe_reg, text='Errors will occur in case of high coupling coefficient and frequency detuning',foreground = 'red')
        # self.warninglabel.place(x=10, y=770)
        self.calclabel = ttk.Label( text='Calculation scheme (non stochastic)')
        self.calclabel.place(x = 760, y = 750)

        self.calctype = ('Euler', 'Runge_Kutta')
        self.calc_var = tk.StringVar()
        self.calc_menu = ttk.OptionMenu(
            self.parent,
            self.calc_var,
            self.calctype[0],
            *self.calctype,
            command=self.parent.filter_changed)
        self.calc_menu.place(x = 760, y = 780)
        
        
        self.Euler_mod = Optionbox(self.parent,('Regular', 'Fast'))
        self.Euler_mod.option_var.set("Regular")
        #self.Euler_mod.bind("<MouseWheel>",  self.empty_scroll_command) # prevent mouse scrolling
        #self.Euler_mod.calc_menu.place ( x = 760, y = 700)
        self.Euler_mod.calc_menu.configure(state = 'disabled')
        self.Euler_mod.calc_menu.place_forget()
        
        self.injfix = Checkbox(self.parent, "Injection multiplication 'fix' ")    
        self.injfix.checkboxm.place(x = 760, y = 810 + 90)
        self.injfix.checkboxm.configure(state='disable')
        self.injfix.checkboxm.place_forget()
        
        self.scipym = Checkbox(self.parent, "Use SCIPY RK45 ") 
        self.scipym.var.set(1)
        self.scipym.checkboxm.place(x = 760, y = 780 + 50)
        self.scipym.checkboxm.configure(state='disable')
        self.scipym.checkboxm.place_forget()
        #self.button = ttk.Button(self.pulseframe_reg,text='clear plot', command=self.parent.clearme, cursor="hand2")
        #self.button.place(x=0, y=800, width=150, height=50)

        #self.filename = tk.Entry()
        #self.filename.place(x = 10, y = 730)

        loadbutton_x = 10
        loadbutton_y = 720
        self.load = ttk.Button(text='Download data',
                               command=self.parent.download, cursor="hand2")
        self.load.place(x=loadbutton_x, y=loadbutton_y)

        self.configb = ttk.Button(
            text='Save parameters', command=self.parent.config_b, cursor="hand2")
        self.configb.place(x=loadbutton_x, y=loadbutton_y-50)

        self.loadconfigb = ttk.Button(
            text='Load parameters', command=self.parent.config_setup, cursor="hand2")
        self.loadconfigb.place(x=loadbutton_x, y=loadbutton_y - 100)
        self.option_var = tk.StringVar()

        self.s1 = tk.BooleanVar()
        self.s1.set(0)
        self.s_plt = ttk.Checkbutton(
            text="Save as csv file", variable=self.s1, onvalue=1, offvalue=0)
        self.s_plt.place(x=loadbutton_x, y=loadbutton_y + 60)

        # option menu
        self.options = ('Driver-Master', 'Driver-Master-WDM', 'Driver-Master-Interferometer', 'Driver-Master-WDM-Interferometer',
                        'Driver-Master-Slave', 'Driver-Master-Slave-Interferometer', 'Driver-Master-Slave-WDM',
                        'Driver-Master-Slave-WDM-Interferometer', 'Driver-Master-WDM-Slave-WDM', 'Driver-Master-WDM-Slave-WDM-Interferometer')
        menu_x = 230
        menu_y = 40
        self.option_menu = ttk.OptionMenu(
            self.parent,
            self.option_var,
            self.options[0],
            *self.options,
            command=self.parent.option_changed)

        self.option_menu.place(x=menu_x, y=menu_y)
        self.label = ttk.Label(text='Select setup:')
        self.label.place(x=menu_x, y=menu_y - 25)

        self.Limage0 = ImageTk.PhotoImage(
            Image.open("config/graphics\L2.png").resize((200, 100)))
        self.parent.Imlbl[0] = ttk.Label(image=self.Limage0)

        self.Limage1 = ImageTk.PhotoImage(
            Image.open("config/graphics\LW.png").resize((220, 100)))
        self.parent.Imlbl[1] = ttk.Label(image=self.Limage1)

        self.Limage2 = ImageTk.PhotoImage(
            Image.open("config/graphics\LI.png").resize((260, 100)))
        self.parent.Imlbl[2] = ttk.Label(image=self.Limage2)

        self.Limage3 = ImageTk.PhotoImage(
            Image.open("config/graphics\LWI.png").resize((300, 100)))
        self.parent.Imlbl[3] = ttk.Label(image=self.Limage3)

        self.Limage4 = ImageTk.PhotoImage(
            Image.open("config/graphics\MS.png").resize((200, 140)))
        self.parent.Imlbl[4] = ttk.Label(image=self.Limage4)

        self.Limage5 = ImageTk.PhotoImage(
            Image.open("config/graphics\MSI.png").resize((260, 140)))
        self.parent.Imlbl[5] = ttk.Label(image=self.Limage5)

        self.Limage6 = ImageTk.PhotoImage(
            Image.open("config/graphics\MSW.png").resize((220, 130)))
        self.parent.Imlbl[6] = ttk.Label(image=self.Limage6)

        self.Limage7 = ImageTk.PhotoImage(
            Image.open("config/graphics\MSWI.png").resize((300, 140)))
        self.parent.Imlbl[7] = ttk.Label(image=self.Limage7)

        self.Limage8 = ImageTk.PhotoImage(
            Image.open("config/graphics\MWSW.png").resize((230, 135)))
        self.parent.Imlbl[8] = ttk.Label(image=self.Limage8)

        self.Limage9 = ImageTk.PhotoImage(Image.open(
            "config/graphics\MWSWI2.png").resize((300, 140)))
        self.parent.Imlbl[9] = ttk.Label(image=self.Limage9)

        self.notebook = ttk.Notebook(self.parent)
        self.notebook.place(x=self.parent.pltsx, y=self.parent.pltsy)
        # create frames
        self.clc_frame = ttk.Frame(self.notebook, width=215, height=560)
        self.plt_frame = ttk.Frame(self.notebook, width=215, height=560)
        self.clc_frame.place(x=self.parent.pltsx, y=self.parent.pltsy)
        self.plt_frame.place(x=self.parent.pltsx, y=self.parent.pltsy)
        # add frames to notebook
        self.notebook.add(self.clc_frame, text='Calculation setup')
        self.notebook.add(self.plt_frame, text='Plot setup')
        self.file_data = np.loadtxt("config/calc_settings.txt")

        self.graphic_output = Checkbox(self.parent, "Turn on external display")    
        self.graphic_output.checkboxm.place(x=loadbutton_x, y=loadbutton_y + 100)
        
        self.power_output = Checkbox(self.parent, "Display laser output in W")
        self.power_output.checkboxm.place(x = loadbutton_x, y = loadbutton_y + 140)
        #self.graphic_output.checkboxm.configure(state='disable')
        names = ["Driver output","Laser output (N,Q)","Phase Calculation","Frequency Calculation",\
                 "Temperature effects","Noise","Biexponential temperature","Interference","Transfer function"]              
        self.master_clcoptions = 9                 #number of master plot options
        self.m_checkbox_clc = [Checkbox(self.clc_frame, names[i],self.clcchange) for i in range(0,self.master_clcoptions)]  
        c = 10
        self.file_data = np.loadtxt("config/calc_settings.txt")
        for i in range(0, self.master_clcoptions):
            self.m_checkbox_clc[i].var.set(int(self.file_data[i]))
            
        for i in range(0, self.master_clcoptions):
            self.m_checkbox_clc[i].checkboxm.place(x = 10, y = 10 + c)
            c = c + 30
 


        namesplt = ["Driver output","Laser output (N)","Laser output (Q)","Phase Calculation",\
                 "Frequency Calculation","Temperature effects","Interference","Transfer function"]
        self.master_pltoptions = 8
        self.m_checkbox_plt = [Checkbox(self.plt_frame, namesplt[i],self.pltchange) for i in range(0, self.master_pltoptions) ]
        self.file_datap = np.loadtxt("config/plt_settings.txt")
        for i in range(0, self.master_pltoptions):
            self.m_checkbox_plt[i].var.set(int(self.file_datap[i]))  
        c = 10
        for i in range(0, self.master_pltoptions):
            self.m_checkbox_plt[i].checkboxm.place (x = 10, y = 10 + c)
            c = c + 30
        
        self.slave_clcoptions = 8  
        self.s_checkbox_clc = [Checkbox(self.clc_frame, names[i]+" S",self.clcchange) for i in range(0, self.slave_clcoptions)]
        for i in range(self.master_clcoptions, self.slave_clcoptions + self.master_clcoptions):
            self.s_checkbox_clc[i-self.master_clcoptions].var.set(int(self.file_data[i])) 

         
        self.slave_pltoptions = 7 
        self.s_checkbox_plt = [Checkbox(self.plt_frame, namesplt[i] + " S",self.pltchange) for i in range(0,self.slave_pltoptions)]
        for i in range( self.master_pltoptions,  self.master_pltoptions + self.slave_pltoptions):
            self.s_checkbox_plt[i-self.master_pltoptions].var.set(int(self.file_datap[i])) 

       
        labelr = ["Transfer function range",""]
        self.transfer_range = [Input(self.parent,0,200,0.1,0, label = labelr[i], tip = "Modulation frequency range for transfer function")
                               for i in range(0,2)]
        
        if (os.path.exists("config/default_config.json")):  # load data from config file
            f = open("config/default_config.json")
            data = json.load(f) 
            self.transfer_range[0].var.set(data["transfer_settings"]["freq_range_start"])
            self.transfer_range[1].var.set(data["transfer_settings"]["freq_range_end"])
            f.close()     
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, len(self.stabinput)):
                self.stabinput[i].var.set("0")
                
        for i in range(len(self.transfer_range)):
                    self.parent.transdata[i] = self.transfer_range[i].var.get()
                    
        self.transfer_range[0].label.place(x=self.parent.wdmx -30,y=self.parent.pltsy + 344)
        for i in range(0,len(self.transfer_range)):
            self.transfer_range[i].spinboxm.bind("<FocusOut>", self.parent.input_reset)
            #self.transfer_range[i].spinboxm.place(x=self.parent.pltsx + 670,y=self.parent.pltsy + 420 + i*40 + 5)
            self.transfer_range[i].spinboxm.place(x=self.parent.wdmx -30,y=self.parent.pltsy + 367 + i*40 + 5)
          
        
    def clcchange(self):
        calc_arr = [self.m_checkbox_clc[i].var.get() for i in range(0,self.master_clcoptions)]+ \
                    [self.s_checkbox_clc[i].var.get() for i in range(0,self.slave_clcoptions)]
        np.savetxt("config/calc_settings.txt",calc_arr)
    def pltchange(self):
        plt_arr = [self.m_checkbox_plt[i].var.get() for i in range(0,self.master_pltoptions)] + \
                    [self.s_checkbox_plt[i].var.get() for i in range(0,self.slave_pltoptions)]
        np.savetxt("config/plt_settings.txt",plt_arr)
        
    def initModFrame(self):
        arbpulsebutton_x = 0
        scrollbarx = 500
        scrollbary = 920
        modframex = 900
        modframey = 960
        self.pulseframe_mod = tk.Frame(self.pulseframe, width = modframex, height = modframey)
        self.pulseframe.add(self.pulseframe_mod, text='Modified Pulses')
        arbpulsebutton_x = 0
        arbpulsebutton_y = 620

        self.mod_frame_canvas = tk.Canvas(self.pulseframe_mod, width=510, height=960)
        self.mod_frame_canvas.place(x = 510, y = 0)
        
        self.scroll_x = ttk.Scrollbar(self.pulseframe_mod, orient = tk.HORIZONTAL, command = self.mod_frame_canvas.xview)
        self.scroll_x.place(x = scrollbarx , y = scrollbary,width = 400, height =20)
        
        self.mod_frame_canvas.configure(xscrollcommand = self.scroll_x.set)
        self.mod_frame_canvas.bind('<Configure>',lambda e: self.mod_frame_canvas.configure(scrollregion = self.mod_frame_canvas.bbox("all")))
        
        def _on_mouse_wheel(event):
            if (self.parent.pulseid  > 3 or self.parent.slavepulseid > 3):
                #self.mod_frame_canvas.xview_scroll(-1 * int((event.delta/90)), "units")
                n = -event.delta / abs(event.delta)     #1 inversion
                p = self.scroll_x.get()[0] + (n*.01) #return scrollbar position and adjust it by a fraction
                self.mod_frame_canvas.xview_moveto(p)
        
        self.mod_frame_canvas.bind_all ("<MouseWheel>", _on_mouse_wheel)
        
        self.pulseframe_mod_aux = tk.Frame( self.mod_frame_canvas, width= 510, height = modframey)#, bg = 'green')
        self.pulseframe_mod_aux.place(x = 0 , y = 0)
        
        self.mod_frame_canvas.create_window((0,0), window = self.pulseframe_mod_aux, anchor = 'nw')
        
        self.pulsebutton = ttk.Button(
            self.pulseframe_mod, text='Build plot', command=self.parent.arbpulsecalculation)
        self.pulsebutton.place(
            x=arbpulsebutton_x, y=arbpulsebutton_y, width=150, height=50)

        self.addpulsebutton = ttk.Button(
            self.pulseframe_mod, text='Add pulse ', command=self.parent.arbpulseadd)
        self.addpulsebutton.place(x=arbpulsebutton_x, y=arbpulsebutton_y + 160, width=120)

        self.addslavepulsebutton = ttk.Button(
            self.pulseframe_mod, text='Add slave pulse ', command=self.parent.slavearbpulseadd)
        self.addslavepulsebutton.place(
            x=arbpulsebutton_x + 130, y=arbpulsebutton_y + 160, width=160)

        self.clearpulsebutton = ttk.Button(
            self.pulseframe_mod, text='Remove pulse', command=self.parent.arbpulseclear)
        self.clearpulsebutton.place(
            x=arbpulsebutton_x, y=arbpulsebutton_y + 210, width=120)

        self.clearslavepulsebutton = ttk.Button(
            self.pulseframe_mod, text='Remove slave pulse', command=self.parent.slavearbpulseclear)
        self.clearslavepulsebutton.place(
            x=arbpulsebutton_x + 130, y=arbpulsebutton_y + 210, width=160)
        
        self.load_trainbutton = ttk.Button(
            self.pulseframe_mod, text='Load pulse train ', command=self.parent.loadconfig_train)
        self.load_trainbutton.place(x=arbpulsebutton_x, y=arbpulsebutton_y + 60)
        
        self.save_trainbutton = ttk.Button(
            self.pulseframe_mod, text='Save pulse train ', command=self.parent.config_train)
        self.save_trainbutton.place(x=arbpulsebutton_x, y=arbpulsebutton_y + 100)
        
        self.spectrum_on = Checkbox(self.pulseframe_mod,"Calculate spectrum of 1 pulse")    
        self.spectrum_on.checkboxm.place(x=arbpulsebutton_x + 260, y=arbpulsebutton_y + 45)
        
        self.spectrum2_on = Checkbox(self.pulseframe_mod,"Calculate spectrum of 2 pulse")    
        self.spectrum2_on.checkboxm.place(x=arbpulsebutton_x + 260, y=arbpulsebutton_y + 85)
        
        
        self.spectrumall_on = Checkbox(self.pulseframe_mod,"Calculate spectrum of the train")    
        self.spectrumall_on.checkboxm.place(x=arbpulsebutton_x + 260, y=arbpulsebutton_y + 125)
        
        
        self.spectrum_menu = Optionbox(self.pulseframe_mod,('Master', 'Slave') )
        self.spectrum_menu.calc_menu.place (x=arbpulsebutton_x + 300, y=arbpulsebutton_y + 160)
        
        dir_path = r'config/trains'
        res = []
        res.append("New train")
        for path in os.listdir(dir_path):
        # check if current path is a file
            if os.path.isfile(os.path.join(dir_path, path)):
                res.append(path)
        self.train_menu = Optionbox(self.pulseframe_mod,res , change = self.parent.train_change, tip = "Select saved train from \ config \ trains folder\nTo clean or make a new train, select - New train ")
        self.train_menu.calc_menu.place (x=arbpulsebutton_x + 260, y=arbpulsebutton_y, width=200)
        
        
        
        
        self.specpulse1 = Input(self.pulseframe_mod, 0, 20, 1, 1)
        self.specpulse1.spinboxm.config(width = 10)
        self.specpulse1.spinboxm.place(x=arbpulsebutton_x + 390, y=arbpulsebutton_y + 160)
        
        self.specpulse2 = Input(self.pulseframe_mod, 0, 20, 1, 1)
        self.specpulse2.spinboxm.config(width = 10)
        self.specpulse2.spinboxm.place(x=arbpulsebutton_x + 390, y=arbpulsebutton_y + 210)
        
        self.pulseidlbl = ttk.Label(self.pulseframe_mod, text = f"Number of master pulses : {self.parent.pulseid}")
        self.pulseidlbl.place(x=arbpulsebutton_x, y=arbpulsebutton_y + 250)
        
        self.slavepulseidlbl =  ttk.Label(self.pulseframe_mod, text = f"Number of slave pulses : {self.parent.slavepulseid}")
        self.slavepulseidlbl.place(x=arbpulsebutton_x, y=arbpulsebutton_y + 290)
        
    def initStabFrame(self):
        values = ['0',str(8*10**3),'-0.5','3.14','18000','100']
        values_start = [0,0,-100,-100,-100,0]
        values_end = [100000,100000,100,100,100000,500]
        increments = [1, 1, 0.01, 0.01, 0.1,1]
        numberofw = 7        
        self.pulseframe_stab = ttk.Frame(self.pulseframe, width = 900, height = 900)
        self.pulseframe.add(self.pulseframe_stab, text='Stability diagram')
        names =['Starting Qinj','Max Qinj','Min phi','Max phi','Find soultion near','Number of points']
        tips =["Initial value of Q", "Maximum value of Q", "Minimum phase value", "Maximum phase value", " Approximate value to find the solution around", "Number of points in the diagram "]
        images = ["Q_min.png","Q_max.png","phi_min.png","phi_max.png","Solve_near.png", "N_points.png"]
        path = r'config/graphics/latex/'
        self.stabinput = [Input(self.pulseframe_stab, values_start[i], values_end[i], increments[i], values[i], pic = path + images[i], tip = tips[i]) for i in range(6)]
        
        self.stabbutton = ttk.Button(self.pulseframe_stab, text = 'Build diagram',command = self.parent.stabcalc, cursor = "hand2")
        self.stabbutton.place(x = 0, y = 800, width = 150, height = 50)
        
        if (os.path.exists("config/default_config.json")):  # load data from config file
            f = open("config/default_config.json")
            data = json.load(f) 
            self.stabinput[0].var.set(data["diagram_settings"]["Qinj_min"])
            self.stabinput[1].var.set(data["diagram_settings"]["Qinj_max"])
            self.stabinput[2].var.set(data["diagram_settings"]["phi_min"])
            self.stabinput[3].var.set(data["diagram_settings"]["phi_max"])
            self.stabinput[4].var.set(data["diagram_settings"]["sol_point"])
            self.stabinput[5].var.set(data["diagram_settings"]["pointsnum"])
            f.close()     
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, len(self.stabinput)):
                self.stabinput[i].var.set("0")
        cc = 0
        for i in range(0,6):
            self.stabinput[i].label.place(x = 400, y = 580 + cc)
            self.stabinput[i].spinboxm.place(x = 500, y = 580 + cc )
            cc = cc + 35
        for i in range(0,len(values)):
            self.parent.stabdata[i] = self.stabinput[i].var.get()
            
    def initDriver(self):
        self.fl = ttk.Label(text='Driver filter type:')
        self.fl.place(x=self.parent.pltsx + 220, y=self.parent.pltsy + 505)

        self.filtertype = ('Bessel', 'Butterworth')
        self.f_var = tk.StringVar()
        self.filter_menu = ttk.OptionMenu(
            self.parent,
            self.f_var,
            self.filtertype[1],
            *self.filtertype,
            command=self.parent.filter_changed)
        self.filter_menu.place(x=self.parent.pltsx + 350,
                               y=self.parent.pltsy + 500)

    
        self.siglbl = ttk.Label(text='Driver signal type:')
        self.siglbl.place(x=self.parent.pltsx + 220, y=self.parent.pltsy + 455)
        
        self.sigtype_menu = Optionbox(self.parent,('Exp', 'Cos') )
        self.sigtype_menu.calc_menu.place (x=self.parent.pltsx + 350, y=self.parent.pltsy + 450)
        
        
        self.drnotebook = ttk.Notebook(self.parent)
        self.drnotebook.place(x=self.parent.pltsx + 220,
                              y=self.parent.pltsy + 70)
        self.maindrframe = ttk.Frame(self.drnotebook, width=240, height = 340)
        self.slavedrframe = ttk.Frame(self.drnotebook, width=240, height = 340)
        self.drnotebook.add(self.maindrframe, text='Main Driver')
        self.drnotebook.add(self.slavedrframe, text='Slave Driver')
        
        if (os.path.exists("config/default_config.json")):
            f = open("config/default_config.json")
            data = json.load(f)
            self.parent.myvar[0].set(data["driver_settings"]["Ibias"]*1e12)
            
            self.parent.myvar[1].set(data["driver_settings"]["Iptp"]*1e12)
            self.parent.myvar[2].set(
                data["driver_settings"]["number of Pulses"])
            self.parent.myvar[3].set(
                data["driver_settings"]["integration rate"])
            self.parent.myvar[4].set(
                data["driver_settings"]["filter bandwidth"])
            self.parent.myvar[5].set(data["driver_settings"]["filter order"])
            self.parent.myvar[6].set(
                data["driver_settings"]["pulse repetiotion rate"])
            self.parent.myvar[7].set(data["driver_settings"]["pulse width"])
            self.parent.myvar[8].set(data["driver_settings"]["decay time"])
            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, self.parent.datadend):
                self.parent.myvar[i].set("0")
        step = [0.1, 0.1, 1, 0.00001, 1, 1, 0.1, 0.1, 0.01]
        limit = [1000,1000,100,1,20,20,20,10,10]
        labels = ["Bias Current (in mA)","PeaktoPeak Current (in mA)","Number of Pulses","Integration rate","Filter bandwidth",
                  "Filter order","Driver pulse repetiotion rate, GHz","The width of the pump current pulse, ns",
                  "Decay time of the electrical circuit, ns"]
        tips = ["Bias Current (in mA)","PeaktoPeak Current (in mA)","Number of Pulses","Integration rate \nLower integration rate makes calculation longer","Filter bandwidth",
                  "Filter order","Driver pulse repetiotion rate, GHz","The width of the pump current pulse, ns \nMaximum width of the pulse is 1/frequency",
                  "Decay time of the electrical circuit, ns"]
        images = ["I_bias.png","I_ptp.png","N_pulses.png","int_rate.png","band.png", "order.png","rep.png","delta_tau.png","tau_d.png"]
        path = r'config/graphics/latex/'
        self.master_entry = [Input(self.maindrframe,0,limit[i],step[i],self.parent.myvar[i].get(),label = "", tip = tips[i],pic = path + images[i] )\
                            for i in range(0,self.parent.datalstart)]
        self.testlbl = ttk.Label(self.maindrframe)
        file = "config/graphics/latex/alpha.png"
        # for i in range(0, self.parent.datadend):
        #       #img =  tk.PhotoImage(file=path + images[i])
        #       self.master_entry[i].label.config(image = self.master_entry[i].image ,text = "", compound = 'left')
        #       self.master_entry[i].label.update()  
        for i in range(0, self.parent.datadend):
            if (self.master_entry[i].var.get() != ''):
                self.parent.data[i] = self.master_entry[i].var.get()
                
        c = 10
        for i in range(0, self.parent.datadend):
            self.master_entry[i].label.place(x = 5, y = c + 3)
            self.master_entry[i].spinboxm.place(x = 120, y= c, width=100)
            c = c + 35
        # c = 20
        # for i in range(0, self.parent.datadend):
        #     self.driverentry[i].spinboxm.place(x=0, y=c, width=100)
        #     # self.parent.entry[i].place(x=0, y=c, width=100)
        #     c = c + 47

    def initSlaveDriver(self):
        if (os.path.exists("config/default_config.json")):
            f = open("config/default_config.json")
            data = json.load(f)
            self.parent.slavevar[0].set(
                data["slave_driver_settings"]["Ibias"]*1e12)
            self.parent.slavevar[1].set(
                data["slave_driver_settings"]["Iptp"]*1e12)
            self.parent.slavevar[2].set(
                data["slave_driver_settings"]["number of Pulses"])
            self.parent.slavevar[3].set(
                data["slave_driver_settings"]["integration rate"])
            self.parent.slavevar[4].set(
                data["slave_driver_settings"]["filter bandwidth"])
            self.parent.slavevar[5].set(
                data["slave_driver_settings"]["filter order"])
            self.parent.slavevar[6].set(
                data["slave_driver_settings"]["pulse repetiotion rate"])
            self.parent.slavevar[7].set(
                data["slave_driver_settings"]["pulse width"])
            self.parent.slavevar[8].set(
                data["slave_driver_settings"]["decay time"])
            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, self.parent.datadend):
                self.parent.myvar[i].set("0")
        
        step = [0.1, 0.1, 1, 0.00001, 1, 1, 0.1, 0.1, 0.01]
        limit = [1000,1000,100,1,20,20,20,10,10]
        labels = ["Bias Current (in mA)","PeaktoPeak Current (in mA)","Number of Pulses","Integration rate  \nLower integration rate makes calculation longer","Filter bandwidth",
                  "Filter order","Driver pulse repetiotion rate, GHz","The width of the pump current pulse, ns",
                  "Decay time of the electrical circuit, ns"]
        tips = ["Bias Current (in mA)","PeaktoPeak Current (in mA)","Number of Pulses","Integration rate","Filter bandwidth",
                  "Filter order","Driver pulse repetiotion rate, GHz","The width of the pump current pulse, ns",
                  "Decay time of the electrical circuit, ns"]
        images = ["I_bias.png","I_ptp.png","N_pulses.png","int_rate.png","band.png", "order.png","rep.png","delta_tau.png","tau_d.png"]
        path = r'config/graphics/latex/'
        
        self.slave_entry = [Input(self.slavedrframe,0,limit[i],step[i],self.parent.slavevar[i].get(), label = "", tip = tips[i], pic = path + images[i])\
                            for i in range(0,self.parent.datalstart)]
        self.slave_entry[3].spinboxm.configure(state = "disabled")
        self.slave_entry[3].spinboxm.bind("<MouseWheel>",  self.empty_scroll_command)
        for i in range(0, self.parent.datadend):
            if (self.slave_entry[i].var.get() != ''):
                self.parent.slavedata[i] = self.slave_entry[i].var.get()

        
    def initLaser(self):
        lasnote_x = 230
        lasnote_y = 550
        self.lasnotebook = ttk.Notebook(self.parent)
        self.lasnotebook.place(x=lasnote_x, y=lasnote_y)
        
        self.mainlasframe = ttk.Frame(self.lasnotebook, width=520, height=386)
        self.slavelasframe = ttk.Frame(self.lasnotebook, width=520, height=386)
        self.lasnotebook.add(self.mainlasframe, text='Main laser')
        self.lasnotebook.add(self.slavelasframe, text='Slave laser')
        self.lnotebook = ttk.Notebook(self.mainlasframe)
        lasnote_in_x = 3
        lasnote_in_y = 3
        self.lnotebook.place(x=lasnote_in_x, y=lasnote_in_y)
        self.lnotebooks = ttk.Notebook(self.slavelasframe)
        self.lnotebooks.place(x=lasnote_in_x, y=lasnote_in_y)
        # create frames
        self.framel1 = ttk.Frame(self.lnotebook, width=500, height=360)
        self.framel2 = ttk.Frame(self.lnotebook, width=500, height=360)
        self.framel1.place(x=0, y=0)
        self.framel2.place(x=0, y=0)
        # add frames to notebook
        self.lnotebook.add(self.framel1, text='General Parameters')
        self.lnotebook.add(self.framel2, text='Thermal and injection')

        self.frames1 = ttk.Frame(self.lnotebooks, width=500, height=360)
        self.frames2 = ttk.Frame(self.lnotebooks, width=500, height=360)
        self.frames1.place(x=0, y=0)
        self.frames2.place(x=0, y=0)
        self.lnotebooks.add(self.frames1, text='General Parameters')
        self.lnotebooks.add(self.frames2, text='Thermal and injection')
        #threshold current data on slave and master lasers
        self.I_th_imag = tk.PhotoImage(file='config/graphics/latex/I_th.png')
        self.threscurrM = Input(self.framel1,0,0,0,0,label = "", tip = "Threshold current value of the master laser")
        self.threscurrS = Input(self.frames1,0,0,0,0,label = "", tip = "Threshold current value of the slave laser")
        self.threscurrM.spinboxm.place(x = 340 , y = 10)
        self.threscurrS.spinboxm.place(x = 340 , y = 10)
        
        self.threscurrM.label.config(image = self.I_th_imag, compound = 'left')
        self.threscurrS.label.config(image = self.I_th_imag, compound = 'left')
        self.threscurrM.label.place(x = 250 , y = 10)
        self.threscurrS.label.place(x = 250 , y = 10)
        
        self.threscurrM.spinboxm.config(state = "disabled")
        self.threscurrS.spinboxm.config(state = "disabled")

        self.threscurrM.spinboxm.bind("<MouseWheel>",  self.empty_scroll_command) # prevent mouse scrolling
        self.threscurrS.spinboxm.bind("<MouseWheel>",  self.empty_scroll_command) # prevent mouse scrolling
        
        if (os.path.exists("config/default_config.json")):  # load data from config file
            f = open("config/default_config.json")
            data = json.load(f)
            self.parent.myvar[9].set(data["laser_settings"]["eta"])
            self.parent.myvar[10].set(data["laser_settings"]["t_ph"])
            self.parent.myvar[11].set(data["laser_settings"]["N_transparency"])
            self.parent.myvar[12].set(data["laser_settings"]["N_threshhold"])
            self.parent.myvar[13].set(data["laser_settings"]["chi"])
            self.parent.myvar[14].set(data["laser_settings"]["C_sp"])
            self.parent.myvar[15].set(data["laser_settings"]["Gamma"])
            self.parent.myvar[16].set(data["laser_settings"]["tau_e"])
            self.parent.myvar[17].set(data["laser_settings"]["alpha"])
            self.parent.myvar[18].set(data["laser_settings"]["rh"])
            self.parent.myvar[19].set(data["laser_settings"]["th"])
            self.parent.myvar[20].set(data["laser_settings"]["lambda"])
            self.parent.myvar[21].set(data["laser_settings"]["kap_om"])
            self.parent.myvar[22].set(data["laser_settings"]["t1"])
            self.parent.myvar[23].set(data["laser_settings"]["w_M"])
            self.parent.myvar[24].set(data["laser_settings"]["om_th"])
            self.parent.myvar[25].set(data["laser_settings"]["Ll"])
            self.parent.myvar[26].set(data["laser_settings"]["kap_inj"])
            self.parent.myvar[27].set(data["laser_settings"]["r1"])
            self.parent.myvar[28].set(data["laser_settings"]["r2"])
            self.parent.myvar[29].set(data["laser_settings"]["rhs"])
            self.parent.myvar[30].set(data["laser_settings"]["k12"])
            self.parent.myvar[31].set(data["laser_settings"]["tauj"])
            self.parent.myvar[32].set(data["laser_settings"]["lambdas"])
            self.parent.myvar[33].set(data["laser_settings"]["Cj"])
            self.parent.myvar[34].set(data["laser_settings"]["int_delay"])
            self.parent.myvar[35].set(data["laser_settings"]["phi1"])
            self.parent.myvar[36].set(data["laser_settings"]["phi2"])
            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(self.parent.datalstart, self.parent.datalend):
                self.parent.myvar[i].set("0")
        tips = ["Differential quantum output", "Photon lifetime inside the cavity, ns","Carrier numbers at transparency ",
                  "Carrier numbers at the threshold", "Gain compression factor, 1/W","Fraction of spontaneously emitted photons that end up in the active mode",
                  "Confinement factor","Effective lifetime of the electron, ns","The linewidth enhancement factor (the Henry factor)",
                  "Thermal resistance between the laser active layer and the heat sink, K/W","Thermal rise-time, ns", 
                  "Wavelength (usually 1,5 micrometers), um","The temperature coefficient of frequency, GHz/K",
                  "The field transmission coefficient of the facet through which \nthe external radiation is injected into the laser cavity" ,
                  "Optical frequency detuning", "omega_th", "Resonator length, micrometeres",
                  "Coupling coefficient \ndo not set too high, may result in errors", "Temperature related phenomenological constant", "Temperature related phenomenological constant" ,"rhs",
                  "Temperature related phenomenological constant", "Temperature related phenomenological constant", "Additional wavelength for temperature", "Heat capacity"]       
        self.eta_imag = tk.PhotoImage(file = "config/graphics/latex/eta.png")        
        self.parent.lbl[9] = ttk.Label(
            self.framel1, image = self.eta_imag, text="")
        self.parent.entry[9] = ttk.Spinbox(
            self.framel1, from_=0, to=1, increment=0.01, textvariable=self.parent.myvar[9])#, command=self.parent.change)
        
        self.tau_ph_imag = tk.PhotoImage(file = "config/graphics/latex/tau_ph.png")
        self.parent.lbl[10] = ttk.Label(
            self.framel1, image = self.tau_ph_imag, text="", compound = 'left')
        self.parent.entry[10] = ttk.Spinbox(
            self.framel1, from_=0, to=1, increment=0.0001, textvariable=self.parent.myvar[10])#, command=self.parent.change)

        self.N_tr_imag = tk.PhotoImage(file = "config/graphics/latex/N_tr.png")
        self.parent.lbl[11] = ttk.Label(self.framel1, image = self.N_tr_imag, text="")
        self.parent.entry[11] = ttk.Spinbox(
            self.framel1, from_=0, to=20, increment=0.1, textvariable = self.parent.myvar[11])#, command=self.parent.change)
        
        self.N_th_imag = tk.PhotoImage(file='config/graphics/latex/N_th.png')
        self.parent.lbl[12] = ttk.Label(self.framel1, image =  self.N_th_imag, text="")
        self.parent.entry[12] = ttk.Spinbox(
            self.framel1, from_=0, to=20, increment=0.1, textvariable = self.parent.myvar[12])#, command=self.parent.change)

        self.chi_imag = tk.PhotoImage(file = "config/graphics/latex/chi.png")
        self.parent.lbl[13] = ttk.Label(
            self.framel1, image = self.chi_imag, text="",compound='left')
        self.parent.entry[13] = ttk.Spinbox(
            self.framel1, from_=0, to=30, increment=0.1, textvariable=self.parent.myvar[13])#, command=self.parent.change)
        
        self.C_sp_imag = tk.PhotoImage(file='config/graphics/latex/C_sp.png')
        self.parent.lbl[14] = ttk.Label(self.framel1, image = self.C_sp_imag, text="")
        self.parent.entry[14] = ttk.Spinbox(
            self.framel1, from_=0, to=30, increment=0.00001, textvariable=self.parent.myvar[14])#, command=self.parent.change)

        self.Gamma_imag = tk.PhotoImage(file='config/graphics/latex/Gamma.png')
        self.parent.lbl[15] = ttk.Label(
            self.framel1, image = self.Gamma_imag, text="" ,compound='left')
        self.parent.entry[15] = ttk.Spinbox(
            self.framel1, from_=0, to=10, increment=0.01, textvariable=self.parent.myvar[15])#, command=self.parent.change)

        self.tau_e_imag = tk.PhotoImage(file = "config/graphics/latex/tau_e.png")
        self.parent.lbl[16] = ttk.Label(
            self.framel1, image = self.tau_e_imag, text="",compound='left')
        self.parent.entry[16] = ttk.Spinbox(
            self.framel1, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[16])#, command=self.parent.change)

        self.alpha_imag = tk.PhotoImage(file = "config/graphics/latex/alpha.png")
        self.parent.lbl[17] = ttk.Label(
            self.framel1, image =  self.alpha_imag, text=":",compound='left')
        self.parent.entry[17] = ttk.Spinbox(
            self.framel1, from_=-100, to=100, increment=0.1, textvariable=self.parent.myvar[17])#, command=self.parent.change)

        self.r_h_imag = tk.PhotoImage(file='config/graphics/latex/r_h.png')
        self.parent.lbl[18] = ttk.Label(self.framel2, image = self.r_h_imag, text="",compound='left')
        self.parent.entry[18] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[18])#, command=self.parent.change)

        self.tau_h_imag = tk.PhotoImage(file='config/graphics/latex/tau_h.png')
        self.parent.lbl[19] = ttk.Label(self.framel2, image = self.tau_h_imag, text="" , compound = 'left')
        self.parent.entry[19] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[19])#, command=self.parent.change)

        self.lambda_imag = tk.PhotoImage(file='config/graphics/latex/lambda.png')
        self.parent.lbl[20] = ttk.Label(
            self.framel2, image = self.lambda_imag, text="",compound='left')
        self.parent.entry[20] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[20])#, command=self.parent.change)
        
        self.kappa_om_imag = tk.PhotoImage(file='config/graphics/latex/kappa_om.png')
        self.parent.lbl[21] = ttk.Label(self.framel2, image = self.kappa_om_imag, text="",compound='left') #ttk.Label(self.framel2, text="\u03BA "+"_\uAB65 /2pi? temperature coefficient of frequency")
        self.parent.entry[21] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.1, textvariable=self.parent.myvar[21])#, command=self.parent.change)

        self.t_1_imag = tk.PhotoImage(file='config/graphics/latex/t_1.png')
        self.parent.lbl[22] = ttk.Label(
            self.framel2, image = self.t_1_imag,text="",compound='left')
        self.parent.entry[22] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[22])#, command=self.parent.config_train)
        self.parent.entry[22].configure(state='disable')
        
        self.om_inj_imag = tk.PhotoImage(file='config/graphics/latex/om_inj.png')
        self.parent.lbl[23] = ttk.Label(self.framel2,image = self.om_inj_imag, text="", compound = 'left') #image = self.inj_imag)#text="\u03C9"+"_inj")
        self.parent.entry[23] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[23])
        self.parent.entry[23].configure(state='disable')
        self.parent.entry[23].bind("<MouseWheel>",  self.empty_scroll_command)
        
        self.Om_th_imag = tk.PhotoImage(file='config/graphics/latex/Om_th.png')
        self.parent.lbl[24] = ttk.Label(self.framel2, image = self.Om_th_imag, text="")
        self.parent.entry[24] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.1, textvariable=self.parent.myvar[24])#, command=self.parent.change)
        self.parent.entry[24].configure(state='disable')
        
        self.parent.lbl[25] = ttk.Label(
            self.framel2, text="L, mm :")
        self.parent.entry[25] = ttk.Spinbox(
            self.framel2, from_=0, to=10, increment=0.0001, textvariable=self.parent.myvar[25])#, command=self.parent.change)
        self.parent.entry[25].configure(state='disable')
        
        self.kappa_inj_imag = tk.PhotoImage(file='config/graphics/latex/kappa_inj.png')
        self.parent.lbl[26] = ttk.Label(self.framel2, image = self.kappa_inj_imag,text = '')#text="\u03BA injection")
        self.parent.entry[26] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.01, textvariable=self.parent.myvar[26])#, command=self.parent.change)
        self.parent.entry[26].configure(state='disable')
        self.parent.entry[26].bind("<MouseWheel>",  self.empty_scroll_command)
        
        self.r_1_imag = tk.PhotoImage(file='config/graphics/latex/r_1.png')
        self.parent.lbl[27] = ttk.Label(self.framel2, image = self.r_1_imag,text = '')
        self.parent.entry[27] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.00001, textvariable=self.parent.myvar[27])#, command=self.parent.change)

        self.r_2_imag = tk.PhotoImage(file='config/graphics/latex/r_2.png')
        self.parent.lbl[28] = ttk.Label(self.framel2, image = self.r_2_imag, text = '')
        self.parent.entry[28] = ttk.Spinbox(
            self.framel2, from_=-10, to=100, increment=0.00001, textvariable=self.parent.myvar[28])#, command=self.parent.change)

        self.r_hs_imag = tk.PhotoImage(file = 'config/graphics/latex/r_hs.png')
        self.parent.lbl[29] = ttk.Label(self.framel2, text="---")
        self.parent.entry[29] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.001, textvariable=self.parent.myvar[29])#, command=self.parent.change)
        self.parent.entry[29].configure(state='disable')
        
        self.k_12_imag = tk.PhotoImage(file = 'config/graphics/latex/k_12.png')
        self.parent.lbl[30] = ttk.Label(self.framel2, image = self.k_12_imag,text = '')
        self.parent.entry[30] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.0001, textvariable=self.parent.myvar[30])#, command=self.parent.change)

        self.tau_j_imag = tk.PhotoImage(file = 'config/graphics/latex/tau_j.png')
        self.parent.lbl[31] = ttk.Label(self.framel2, image = self.tau_j_imag,text = '')
        self.parent.entry[31] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.0001, textvariable=self.parent.myvar[31])#, command=self.parent.change)

        self.parent.lbl[32] = ttk.Label(
            self.framel2,text="---", compound = 'left')
        self.parent.entry[32] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.01, textvariable=self.parent.myvar[32])#, command=self.parent.change)
        self.parent.entry[32].configure(state='disable')
        
        self.C_j_imag = tk.PhotoImage(file = 'config/graphics/latex/C_j.png')
        self.parent.lbl[33] = ttk.Label(self.framel2, image = self.C_j_imag,text = '', compound = "left")
        self.parent.entry[33] = ttk.Spinbox(
            self.framel2, from_=0, to=100, increment=0.01, textvariable=self.parent.myvar[33])#, command=self.parent.change)

        self.mastertip = [Hovertip(self.parent.entry[i], tips[i-self.parent.datalstart], hover_delay=1000)\
                            for i in range(self.parent.datadend,34)]
        self.masterlbltip = [Hovertip(self.parent.lbl[i], tips[i-self.parent.datalstart], hover_delay=1000)\
                            for i in range(self.parent.datadend,34)]
            
        self.unnecessary = [22,24,25,29,32]#hide irrelevant parameteres entries    
        
        for i in range(self.parent.datalstart, self.parent.datalend):
            if (self.parent.entry[i].get() != ''):
                self.parent.data[i] = self.parent.entry[i].get()
        c = 10
        for i in range(9, 18):
            self.parent.entry[i].place(x=90, y=c, width=110)
            self.parent.lbl[i].place(x=0, y=c+3)
            c = c + 35
        c = 10
        
        for i  in range(18, 26):
            if( (i in self.unnecessary) != 1):
                self.parent.entry[i].place(x=150, y=c, width=110)
                self.parent.lbl[i].place(x=0, y=c)
                c = c + 35
        c = 10

        for i in range(26, 34):
            if( (i in self.unnecessary) != 1 ):
                self.parent.entry[i].place(x=370, y=c, width=110)
                self.parent.lbl[i].place(x=320, y=c+3)
                c = c + 35

        self.savebutton = ttk.Button(
            text='Save', command=self.parent.savelaser)
        self.savebutton.place(x=self.parent.pltsx + 750,
                              y=self.parent.pltsy + 630, height=50, width=100)

    def initSlaveLaser(self):
        if (os.path.exists("config/default_config.json")):
            f = open("config/default_config.json")
            data = json.load(f)

            self.parent.slavevar[9].set(data["slave_laser_settings"]["eta"])
            self.parent.slavevar[10].set(data["slave_laser_settings"]["t_ph"])
            self.parent.slavevar[11].set(
                data["slave_laser_settings"]["N_transparency"])
            self.parent.slavevar[12].set(
                data["slave_laser_settings"]["N_threshhold"])
            self.parent.slavevar[13].set(data["slave_laser_settings"]["chi"])
            self.parent.slavevar[14].set(data["slave_laser_settings"]["C_sp"])
            self.parent.slavevar[15].set(data["slave_laser_settings"]["Gamma"])
            self.parent.slavevar[16].set(data["slave_laser_settings"]["tau_e"])
            self.parent.slavevar[17].set(data["slave_laser_settings"]["alpha"])
            self.parent.slavevar[18].set(data["slave_laser_settings"]["rh"])
            self.parent.slavevar[19].set(data["slave_laser_settings"]["th"])
            self.parent.slavevar[20].set(
                data["slave_laser_settings"]["lambda"])
            self.parent.slavevar[21].set(
                data["slave_laser_settings"]["kap_om"])
            self.parent.slavevar[22].set(data["slave_laser_settings"]["t1"])
            self.parent.slavevar[23].set(data["slave_laser_settings"]["w_M"])
            self.parent.slavevar[24].set(data["slave_laser_settings"]["om_th"])
            self.parent.slavevar[25].set(data["slave_laser_settings"]["Ll"])
            self.parent.slavevar[26].set(
                data["slave_laser_settings"]["kap_inj"])
            self.parent.slavevar[27].set(data["slave_laser_settings"]["r1"])
            self.parent.slavevar[28].set(data["slave_laser_settings"]["r2"])
            self.parent.slavevar[29].set(data["slave_laser_settings"]["rhs"])
            self.parent.slavevar[30].set(data["slave_laser_settings"]["k12"])
            self.parent.slavevar[31].set(data["slave_laser_settings"]["tauj"])
            self.parent.slavevar[32].set(
                data["slave_laser_settings"]["lambdas"])
            self.parent.slavevar[33].set(data["slave_laser_settings"]["Cj"])
            self.parent.slavevar[34].set(
                data["slave_laser_settings"]["int_delay"])
            self.parent.myvar[35].set(data["laser_settings"]["phi1"])
            self.parent.myvar[36].set(data["laser_settings"]["phi2"])
            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(self.parent.datalstart, self.parent.datalend):
                self.parent.slavevar[i].set("0")
                
        self.parent.slavelbl[9] = ttk.Label(
            self.frames1, image = self.eta_imag, text=" ", compound = 'left')
        self.parent.slaveentry[9] = ttk.Spinbox(
            self.frames1, from_=0, to=1, increment=0.01, textvariable=self.parent.slavevar[9])#, command=self.parent.slchange)

        self.parent.slavelbl[10] = ttk.Label(
            self.frames1, image = self.tau_ph_imag, text="", compound = 'left')
        self.parent.slaveentry[10] = ttk.Spinbox(
            self.frames1, from_=0, to=1, increment=0.0001, textvariable=self.parent.slavevar[10])#, command=self.parent.slchange)

        self.parent.slavelbl[11] = ttk.Label(
            self.frames1, image = self.N_tr_imag, text=" ", compound = 'left')
        self.parent.slaveentry[11] = ttk.Spinbox(
            self.frames1, from_=0, to=20, increment=0.1, textvariable=self.parent.slavevar[11])#, command=self.parent.slchange)

        self.parent.slavelbl[12] = ttk.Label(
            self.frames1, image = self.N_th_imag, text="" , compound = 'left')
        self.parent.slaveentry[12] = ttk.Spinbox(
            self.frames1, from_=0, to=20, increment=0.1, textvariable=self.parent.slavevar[12])#, command=self.parent.slchange)

        self.parent.slavelbl[13] = ttk.Label(
            self.frames1, image = self.chi_imag, text="",compound = 'left')
        self.parent.slaveentry[13] = ttk.Spinbox(
            self.frames1, from_=0, to=30, increment=0.1, textvariable=self.parent.slavevar[13])#, command=self.parent.slchange)

        self.parent.slavelbl[14] = ttk.Label(
            self.frames1,image = self.C_sp_imag, text="")
        self.parent.slaveentry[14] = ttk.Spinbox(
            self.frames1, from_=0, to=30, increment=0.00001, textvariable=self.parent.slavevar[14])#, command=self.parent.slchange)

        self.parent.slavelbl[15] = ttk.Label(
            self.frames1, text="Г")
        self.parent.slaveentry[15] = ttk.Spinbox(
            self.frames1, from_=0, to=10, increment=0.01, textvariable=self.parent.slavevar[15])#, command=self.parent.slchange)

        self.parent.slavelbl[16] = ttk.Label(
            self.frames1, image = self.tau_e_imag, text="", compound = 'left')
        self.parent.slaveentry[16] = ttk.Spinbox(
            self.frames1, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[16])#, command=self.parent.slchange)

        self.parent.slavelbl[17] = ttk.Label(
            self.frames1, image = self.alpha_imag, text="", compound='left')
        self.parent.slaveentry[17] = ttk.Spinbox(
            self.frames1, from_=-100, to=100, increment=0.1, textvariable=self.parent.slavevar[17])#, command=self.parent.slchange)

        self.parent.slavelbl[18] = ttk.Label(self.frames2,imag = self.r_h_imag, text="", compound='left')
        self.parent.slaveentry[18] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[18])#, command=self.parent.slchange)

        self.parent.slavelbl[19] = ttk.Label(
            self.frames2, image = self.tau_h_imag, text="", compound='left')
        self.parent.slaveentry[19] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[19])#, command=self.parent.slchange)

        self.parent.slavelbl[20] = ttk.Label(
            self.frames2,image = self.lambda_imag, text="",compound='left')
        self.parent.slaveentry[20] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[20])#, command=self.parent.slchange)

        self.parent.slavelbl[21] = ttk.Label(
            self.frames2, image = self.kappa_om_imag, text="", compound = 'left')
        self.parent.slaveentry[21] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.1, textvariable=self.parent.slavevar[21])#, command=self.parent.slchange)

        self.parent.slavelbl[22] = ttk.Label(
            self.frames2, image = self.t_1_imag, text=" ",compound = 'left')
        self.parent.slaveentry[22] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[22])#, command=self.parent.slchange)
        self.parent.slaveentry[22].configure(state='disable')
        
        self.parent.slavelbl[23] = ttk.Label(self.frames2, image = self.om_inj_imag, text="  " , compound = 'left')
        self.parent.slaveentry[23] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[23])#, command=self.parent.slchange)

        self.parent.slavelbl[24] = ttk.Label(self.frames2, image = self.Om_th_imag, text="")
        self.parent.slaveentry[24] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.1, textvariable=self.parent.slavevar[24])#, command=self.parent.slchange)
        self.parent.slaveentry[24].configure(state='disable')
        
        self.parent.slavelbl[25] = ttk.Label(self.frames2, text="Ll, \u03BC" +"m" )
        self.parent.slaveentry[25] = ttk.Spinbox(
            self.frames2, from_=0, to=10, increment=0.0001, textvariable=self.parent.slavevar[25])#, command=self.parent.slchange)
        self.parent.slaveentry[25].configure(state='disable')
        
        self.parent.slavelbl[26] = ttk.Label(self.frames2, image = self.kappa_inj_imag)#ttk.Label(self.frames2, text="\u03BA"+"_inj")
        
        self.parent.slaveentry[26] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.01, textvariable=self.parent.slavevar[26])#, command=self.parent.slchange)

        self.parent.slavelbl[27] = ttk.Label(self.frames2, image = self.r_1_imag, text="")
        self.parent.slaveentry[27] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.00001, textvariable=self.parent.slavevar[27])#, command=self.parent.slchange)

        self.parent.slavelbl[28] = ttk.Label(self.frames2, image = self.r_2_imag, text="")
        self.parent.slaveentry[28] = ttk.Spinbox(
            self.frames2, from_=-10, to=100, increment=0.00001, textvariable=self.parent.slavevar[28])#, command=self.parent.slchange)

        self.parent.slavelbl[29] = ttk.Label(self.frames2, text="---")
        self.parent.slaveentry[29] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.001, textvariable=self.parent.slavevar[29])#, command=self.parent.slchange)
        self.parent.slaveentry[29].configure(state='disable')
        
        self.parent.slavelbl[30] = ttk.Label(self.frames2, image = self.k_12_imag, text="")
        self.parent.slaveentry[30] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.0001, textvariable=self.parent.slavevar[30])#, command=self.parent.slchange)

        self.parent.slavelbl[31] = ttk.Label(self.frames2, image = self.tau_j_imag, text=" ")
        self.parent.slaveentry[31] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.0001, textvariable=self.parent.slavevar[31])#, command=self.parent.slchange)

        self.parent.slavelbl[32] = ttk.Label(self.frames2, text="---" , compound = 'left')
        self.parent.slaveentry[32] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.01, textvariable=self.parent.slavevar[32])#, command=self.parent.slchange)
        self.parent.slaveentry[32].configure(state='disable')
        
        self.parent.slavelbl[33] = ttk.Label(self.frames2, image = self.C_j_imag, text="")
        self.parent.slaveentry[33] = ttk.Spinbox(
            self.frames2, from_=0, to=100, increment=0.01, textvariable=self.parent.slavevar[33])#, command=self.parent.slchange)
        
        tips = ["Differential quantum output", "Photon lifetime inside the cavity, ns","Carrier numbers at transparency ",
                  "Carrier numbers at the threshold", "Gain compression factor (in 1/W)","Fraction of spontaneously emitted photons that end up in the active mode",
                  "Confinement factor","Effective lifetime of the electron, ns","The linewidth enhancement factor (the Henry factor)",
                  "Thermal resistance between the laser active layer and the heat sink, K/W","Thermal rise-time, ns", 
                  "Wavelength (usually 1,5 micrometers), um","The temperature coefficient of frequency, GHz/K",
                  "The field transmission coefficient of the facet through which \nthe external radiation is injected into the laser cavity" ,
                  "Optical frequency detuning", "omega_th", "resonator length, micrometeres",
                  "Coupling coefficient \ndo not set too high, may result in errors", "r1", "r2","rhs",
                  "k12", "Tau_j", "Additional wavelength for temperature", "Heat capacity"] 
        
        self.slavetip = [Hovertip(self.parent.slaveentry[i], tips[i-self.parent.datalstart], hover_delay=1000)\
                            for i in range(9,34)]
            
        for i in range(self.parent.datalstart, self.parent.datalend):
            if (self.parent.slaveentry[i].get() != ''):
                self.parent.slavedata[i] = self.parent.slaveentry[i].get()
        #debug - slave entries placement
        # c = 10
        # for i in range(9, 18):
        #     self.parent.slaveentry[i].place(x=0, y=c, width=110)
        #     self.parent.slavelbl[i].place(x=120, y=c)
        #     c = c + 35
        # c = 10

        # for i in range(18, 26):
        #     self.parent.slaveentry[i].place(x=0, y=c, width=110)
        #     self.parent.slavelbl[i].place(x=120, y=c)
        #     c = c + 35
        # c = 10

        # for i in range(26, 34):
        #     self.parent.slaveentry[i].place(x=250, y=c, width=110)
        #     self.parent.slavelbl[i].place(x=370, y=c)
        #     c = c + 35
        # c = 0
        self.slavetip = [Hovertip(self.parent.slaveentry[i], tips[i-self.parent.datalstart], hover_delay=1000)\
                            for i in range(self.parent.datadend,34)]
        self.slavetip = [Hovertip(self.parent.slavelbl[i], tips[i-self.parent.datalstart], hover_delay=1000)\
                            for i in range(self.parent.datadend,34)]
            

    def initWDM(self):
        if (os.path.exists("config/default_config.json")):
            f = open("config/default_config.json")
            data = json.load(f)

            self.parent.wdmvar[0].set(data["wdm_settings"]["filter bandwidth"])
            self.parent.wdmvar[1].set(data["wdm_settings"]["filter order"])
            self.parent.wdmvar[2].set(data["wdm_settings"]["filter shift"])

            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, 3):
                self.parent.wdmvar[i].set("0")
        self.wdmtype = ttk.Label(text='WDM filter type:')
        self.wdmtype.place(x=self.parent.pltsx + 490,
                           y=self.parent.pltsy + 475)

        self.wdmvar1 = tk.StringVar()
        self.wdm_menu = ttk.OptionMenu(
            self.parent,
            self.wdmvar1,
            self.filtertype[1],
            *self.filtertype,
            command=self.parent.filter_changed)
        self.wdm_menu.place(x=self.parent.pltsx + 490,
                            y=self.parent.pltsy + 500)

        self.wdm_imag0 = tk.PhotoImage(file = "config/graphics/latex/band.png")
        self.wdm_imag1 = tk.PhotoImage(file = "config/graphics/latex/order.png")
        self.wdm_imag2 = tk.PhotoImage(file = "config/graphics/latex/shift.png")
        
        self.parent.setuplbl[0] = ttk.Label(image = self.wdm_imag0)#text="Filter Bandwidth ")
        self.parent.dentry[0] = ttk.Spinbox(
            from_=0, to=1000, increment=0.1, textvariable=self.parent.wdmvar[0])#, command=self.parent.savelaser)

        self.parent.setuplbl[1] = ttk.Label(image = self.wdm_imag1)#text="Filter order ")
        self.parent.dentry[1] = ttk.Spinbox(
            from_=0, to=1000, increment=1, textvariable=self.parent.wdmvar[1])#, command=self.parent.savelaser)

        self.parent.setuplbl[2] = ttk.Label(image = self.wdm_imag2)#text="Filter shift")
        self.parent.dentry[2] = ttk.Spinbox(
            from_=-100, to=100, increment=0.1, textvariable=self.parent.wdmvar[2])#, command=self.parent.savelaser)
        for i in range(0,3):
            self.parent.wdmdata[i] =  self.parent.dentry[i].get()  
        
    def initSlaveWDM(self):
        if (os.path.exists("config/default_config.json")):
            f = open("config/default_config.json")
            data = json.load(f)

            self.parent.slavewdmwar[0].set(
                data["slave_wdm_settings"]["filter bandwidth"])
            self.parent.slavewdmwar[1].set(
                data["slave_wdm_settings"]["filter order"])
            self.parent.slavewdmwar[2].set(
                data["slave_wdm_settings"]["filter shift"])

            f.close
        else:
            tk_messagebox.showinfo(
                "Error", "configuration file has not been found!")
            for i in range(0, 3):
                self.parent.slavewdmwar[i].set("0")
        self.wdmslavetype = ttk.Label(text='Slave WDM filter type:')
        self.wdmslavetype.place(x=self.parent.pltsx +
                                630, y=self.parent.pltsy + 475)

        self.wdmvar2 = tk.StringVar()
        self.slavewdm_menu = ttk.OptionMenu(
            self.parent,
            self.wdmvar2,
            self.filtertype[1],
            *self.filtertype,
            command=self.parent.filter_changed)

        self.slavewdm_menu.place(
            x=self.parent.pltsx + 630, y=self.parent.pltsy + 500)
        
        self.wdms_imag0 = tk.PhotoImage(file = "config/graphics/latex/band.png")
        self.wdms_imag1 = tk.PhotoImage(file = "config/graphics/latex/order.png")
        self.wdms_imag2 = tk.PhotoImage(file = "config/graphics/latex/shift.png")
        
        self.parent.slavewdmlbl[0] = ttk.Label(image = self.wdms_imag0)
        self.parent.slavedentry[0] = ttk.Spinbox(
            from_=0, to=1000, increment=0.1, textvariable=self.parent.slavewdmwar[0])

        self.parent.slavewdmlbl[1] = ttk.Label(image = self.wdms_imag1)
        self.parent.slavedentry[1] = ttk.Spinbox(
            from_=0, to=1000, increment=1, textvariable=self.parent.slavewdmwar[1])

        self.parent.slavewdmlbl[2] = ttk.Label(image = self.wdms_imag2)
        self.parent.slavedentry[2] = ttk.Spinbox(
            from_=-100, to=100, increment=0.1, textvariable=self.parent.slavewdmwar[2])
        for i in range(0,3):
            self.parent.slavewdmdata[i] =  self.parent.slavedentry[i].get()
    
    def initInt(self):
        

        self.int_imag0 = tk.PhotoImage(file = "config/graphics/latex/tau_int.png")
        self.int_imag1 = tk.PhotoImage(file = "config/graphics/latex/theta_1.png")
        self.int_imag2 = tk.PhotoImage(file = "config/graphics/latex/theta_2.png")

        self.parent.intlbl[0] = ttk.Label(image = self.int_imag0, text="")

        self.parent.intvar[0].set(self.parent.myvar[34].get())
        self.parent.intentry[0] = ttk.Spinbox(
            from_=0, to=1000, increment=0.1, textvariable=self.parent.intvar[0])#, command=self.parent.savelaser)

        self.parent.intvar[1].set(self.parent.myvar[35].get())
        self.parent.intlbl[1] = ttk.Label(image = self.int_imag1, text="")
        self.parent.intentry[1] = ttk.Spinbox(
            from_=-100, to=1000, increment=1, textvariable=self.parent.intvar[1])#, command=self.parent.savelaser)

        self.parent.intvar[2].set(self.parent.myvar[36].get())
        self.parent.intlbl[2] = ttk.Label(image = self.int_imag2, text="")
        self.parent.intentry[2] = ttk.Spinbox(
            from_=-100, to=100, increment=0.1, textvariable=self.parent.intvar[2])#, command=self.parent.savelaser)
        for i in range(0,3):
            self.parent.intdata[i] =  self.parent.intentry[i].get()
    def InitPhot(self):
        self.parent.photvar[0].set(10)
        self.parent.photentry[0] = ttk.Spinbox(
            from_=0, to=1000, increment=0.1, textvariable=self.parent.photvar[0])
        self.parent.photlbl[0] = ttk.Label(text="Photodetector Bandwidth ")

        self.parent.photvar[1].set(2)
        self.parent.photentry[1] = ttk.Spinbox(
            from_=0, to=1000, increment=1, textvariable=self.parent.photvar[1])
        self.parent.photlbl[1] = ttk.Label(text="Photodetector order ")

        self.parent.photvar[2].set(0)
        self.parent.photentry[2] = ttk.Spinbox(
            from_=0, to=1000, increment=0.1, textvariable=self.parent.photvar[2])
        self.parent.photlbl[2] = ttk.Label(text="Photodetector shift")
        photx = 760
        photy = 440
        c = 0
        for i in range(0, 3):
            self.parent.photentry[i].place(x=self.parent.wdmx+300, y=photy + c + 2, width=100)
            c = c + 48
        c = -40
        for i in range(0, 3):
            self.parent.photlbl[i].place(x=self.parent.wdmx+300, y=photy + 20 + c)
            c = c + 48
        for i in range(0,3):
            self.parent.photdata[i] =  self.parent.photentry[i].get()    