from pandas import DataFrame
import tkinter as tk
from tkinter import messagebox as mb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
from PIL import ImageTk, Image
from idlelib.tooltip import Hovertip

class Input:
    def __init__(self, frame, start , end ,inc , value, change = None, label = None, tip = None, pic = None):
        
        self.var = tk.StringVar()
        self.var.set(value)
        self.spinboxm = ttk.Spinbox( frame, from_ = start, to = end , increment = inc, textvariable= self.var,command = change, width = 10)
        self.spintip = Hovertip(self.spinboxm, tip, hover_delay=1000)
        self.label = ttk.Label(frame, text = label)
        if (pic != None):
            self.image =  tk.PhotoImage(file=pic)
            self.label = ttk.Label(frame, image = self.image, text = label, compound = 'left')
        else:
            self.label = ttk.Label(frame,  text = label)
        self.labeltip = Hovertip(self.label, tip, hover_delay=1000)
        #self.label = ttk.Label(frame, text = label, image = None)
class Checkbox:
    def __init__(self,frame,label , change = None):
        self.var = tk.BooleanVar()
        self.var.set(0)
        self.checkboxm = ttk.Checkbutton(frame, text= label, variable = self.var, onvalue=1, offvalue=0, command = change)
class Optionbox:
    def __init__(self,frame, array, change = None, tip = None):    
        self.options = array
        self.option_var = tk.StringVar()
        self.calc_menu = ttk.OptionMenu(
        frame,
        self.option_var,
        self.options[0],
        *self.options, command = change)
        self.tip = Hovertip(self.calc_menu, tip, hover_delay=1000)         