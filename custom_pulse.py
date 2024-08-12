from pandas import DataFrame
import tkinter as tk
from tkinter import messagebox as mb
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

class Arbpulse:
    def __init__(self,frame = None ):
        self.var = [] #arbitary pulses
        self.data = np.zeros(10) 
        self.lbl =[] 
        self.entry =[]
        self.time = []
        self.filtered = []
        
        for i in range(0,10):
            self.var.append(tk.StringVar())
            self.lbl.append(ttk.Label())
            self.entry.append(ttk.Spinbox())
        self.lbl.append(ttk.Label()) #pulse number
        for i in range(0,10):
                self.var[i].set("0")
        self.var[0].set("6")
        self.lbl[0] = ttk.Label(frame, text="Bias Current, mA")
        self.entry[0] = ttk.Spinbox(frame, from_=-100 , to = 1000, increment= 0.1,textvariable = self.var[0],command = self.change)
            
        self.var[1].set("30")
        self.lbl[1] = ttk.Label(frame, text="PeaktoPeak Current, mA")
        self.entry[1] = ttk.Spinbox(frame, from_=-100 , to = 1000, increment= 0.1,textvariable = self.var[1],command = self.change)
        
        self.var[2].set("0.3")
        self.lbl[2] = ttk.Label(frame, text="Width of the pulse ,ns")
        self.entry[2] = ttk.Spinbox(frame, from_=0 , to = 1000, increment= 0.01,textvariable = self.var[2],command = self.change)

        self.lbl[3] = ttk.Label(frame, text="Gap Current (in mA)")
        self.entry[3] = ttk.Spinbox(frame, from_=-100 , to = 1000, increment= 0.01,textvariable = self.var[3],command = self.change)

        self.lbl[4] = ttk.Label(frame, text="Gap width ns")
        self.entry[4] = ttk.Spinbox(frame, from_=0 , to = 1000, increment= 0.01,textvariable = self.var[4],command = self.change)
        
        self.lbl[5] = ttk.Label(frame, text="Gap position, ns")
        self.entry[5] = ttk.Spinbox(frame, from_=0 , to = 1000, increment= 0.01,textvariable = self.var[5],command = self.change)
        
        self.var[6].set("1.25")
        self.lbl[6] = ttk.Label(frame, text="Repetiton Rate, Ghz")
        self.entry[6] = ttk.Spinbox(frame ,from_=0 , to = 20, increment= 0.01,textvariable = self.var[6], command = self.change)
        
        self.lbl[7] = ttk.Label(frame, text="Time shift, ns")
        self.entry[7] = ttk.Spinbox(frame, from_=0 , to = 20, increment= 0.01,textvariable = self.var[7], command = self.change)
        
        self.var[8].set("0.0001")
        self.lbl[8] = ttk.Label(frame, text="Integration rate")
        self.entry[8] = ttk.Spinbox(frame, from_=0 , to = 20, increment= 0.0001,textvariable = self.var[8], command = self.change)
        
        self.var[9].set("1.00")
        self.lbl[9] = ttk.Label(frame, text="Incline")
        self.entry[9] = ttk.Spinbox(frame, from_=0 , to = 20, increment= 0.01,textvariable = self.var[9], command = self.change)
        
        self.lbl[10] = ttk.Label(frame, text="-")
        
    def put(self,lx,ly,add,pid,frame = None,deltax = None, scroll = None ,canvas = None):
        
        if (pid > 0):
            self.lbl[10].configure(text = f" Pulse {pid}")
            self.lbl[10].place(x = lx, y = ly)
        for i in range(0,10):
            if (self.entry[i].get() != ''):
                self.data[i] = self.entry[i].get()            
        c = ly
        if( pid < 1):
            for i in range(0,10):
                self.lbl[i].place(x =lx, y =  c)
                c = c + 47
        c = ly + 22
        for i in range(0,10):
            self.entry[i].place(x = lx, y =  c , width = 100)
            c = c + 47 
    def hide(self,canvas = None):
        for i in range(0,10):
            self.entry[i].place_forget()
            self.lbl[i].place_forget()
        self.lbl[10].place_forget()
    def change(self):
        for i in range(0,10):
            if (self.entry[i].get() != ''):
                self.data[i] = self.entry[i].get()       
        