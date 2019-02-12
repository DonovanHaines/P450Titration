import tkinter as tk  #to use graphical user interface (windowed program)
#from tkinter import filedialog
#import tkinter.messagebox
import pandas as pd
#from pandastable import Table, TableModel
import pandastable as pdt
import numpy as np
import math
from scipy.optimize import optimize
import lmfit

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import spectra

class TitrationWindow(tk.Toplevel):
    def __init__(self,master=None, spectra=None, reference_spectra=None, width=600, height=800, fit_control_dict=None):
        tk.Toplevel.__init__(self, width=width, height=height)
    self.master=master
    self.loaded_spectra=spectra
    self.reference_spectra=reference_spectra
    
    # Set up the menu system
     menu=tk.Menu(master)
     new_item=tk.Menu(menu)
     new_item.add_command(label="Load Volume Data", command=self.load_volumes)
     new_item.add_command(label="Save Table", command=self.save_table)
     new_item.add_command(label="Exit", command=master.destroy)
     menu.add_cascade(label='File', menu=new_item)
     master.config(menu=menu)
     
     new_item2=tk.Menu(menu)
     new_item2.add_command(label="Simple Fit", command=self.fit_simple)
     new_item2.add_command(label="Cooperative Fit", command=self.fit_coop)

     menu.add_cascade(label="Fit", menu=new_item2)
     
     new_item3=tk.Menu(menu)
     new_item3.add_command(label="About this program", command=self.about)
     menu.add_cascade(label="Help", menu=new_item3)
    
     # Now to add controls
