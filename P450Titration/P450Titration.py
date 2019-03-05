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

import titrationwindow
import spectra



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

class MyGUI(tk.Tk):
    def __init__(self, *args, **kwargs):
      tk.Tk.__init__(self, *args, **kwargs)
      #self.master = master
      #master = self
      self.title("P450 Titration Analysis")
      self.geometry('1920x1080')
      self.loaded_spectra = None
      self.reference_spectra = None
       
      # Set up the menu system
      menu=tk.Menu(self)
      new_item=tk.Menu(menu)
      new_item.add_command(label="Load Spectra", command=self.load_spectra)
      new_item.add_command(label="Load Reference Spectra", command=self.load_reference)
      new_item.add_command(label="Exit", command=self.destroy)
      menu.add_cascade(label='File', menu=new_item)
      self.config(menu=menu)
            
      new_item2=tk.Menu(menu)
      new_item2.add_command(label="Baseline Correct", command=self.baseline_all)
      #new_item2.add_command(label="Baseline Correct All", command=self.baseline_all)
      new_item2.add_command(label="Fit  binding", command=self.fitsimple)
      #new_item2.add_command(label="Fit coop binding", command=self.fitcoop)
      menu.add_cascade(label="Analyze", menu=new_item2)

      new_item3=tk.Menu(menu)
      new_item3.add_command(label="About this program", command=self.about)
      menu.add_cascade(label="Help", menu=new_item3)

          
    ##################################################################################################################################
    ##################################################################################################################################
    def load_spectra(self):
      self.loaded_spectra = loaded_spectra = spectra.spectral_collection(load=True)
      print(loaded_spectra)

      f1=tk.Frame(self, width=400)
      f1.pack(side = tk.LEFT, fill=tk.Y,expand=0)
      df = loaded_spectra.df
      pt = pdt.Table(f1, dataframe=df,
                                    showtoolbar=True, showstatusbar=True)
      pt.show()

      f2 = tk.Frame(self)
      f2.pack(side = tk.TOP, fill = tk.BOTH, expand=1)

      myfig = plt.Figure(figsize=(5,5), dpi=100)
      aplot = myfig.add_subplot(111)
      aplot.plot(loaded_spectra.spectra[0].wavelength,loaded_spectra.spectra[0].absorbance)
      for count in range(1,loaded_spectra.number_spectra):
          aplot.plot(loaded_spectra.spectra[count].wavelength,loaded_spectra.spectra[count].absorbance)
      aplot.set_title("Raw Spectra Loaded")
      canvas = FigureCanvasTkAgg(myfig, f2)
      canvas.draw()
      canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
      toolbar = NavigationToolbar2Tk(canvas, f2)
      toolbar.update()
      canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    ##################################################################################################################################
    def load_reference(self):
        self.reference_spectra = reference_spectra = spectra.spectral_collection(load=True)
        df = reference_spectra.df
      
        f2 = tk.Frame(self)
        f2.pack(side = tk.BOTTOM, fill = tk.BOTH, expand=1)

        myfig = plt.Figure(figsize=(5,5), dpi=100)
        aplot = myfig.add_subplot(111)
        aplot.plot(reference_spectra.spectra[0].wavelength,reference_spectra.spectra[0].absorbance)
        for count in range(1,reference_spectra.number_spectra):
            aplot.plot(reference_spectra.spectra[count].wavelength,reference_spectra.spectra[count].absorbance)
        aplot.set_title("Reference Spectra Loaded")
        canvas = FigureCanvasTkAgg(myfig, f2)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, f2)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    ##################################################################################################################################
    ##################################################################################################################################
    def baseline_all(self):
        fit_all_window = FitAllWindow(master=self, spectra = self.loaded_spectra, reference_spectra = self.reference_spectra)


    ##################################################################################################################################
    ##################################################################################################################################
    def baseline(self):
        result, components, fit_x, fit_y = self.loaded_spectra.spectra[0].baseline_correct(self.reference_spectra, fit_bg=True, fit_scatter=True, fit_gaussian=False, wave_low=300, wave_high=800)
      
        resultwin = tk.Toplevel(width=600, height=800) #make new popup window
        resultwin.wm_title("Result of first fit")
      
        textbuffer = tk.Text(resultwin)
        textbuffer.insert(tk.END, result.fit_report())
        textbuffer.pack(side = tk.TOP, fill = tk.X, expand=1)
        f5 = tk.Frame(resultwin)
        f5.pack(side = tk.BOTTOM, fill = tk.BOTH, expand=1)

        myfig5 = plt.Figure(figsize=(5,5), dpi=100)
        aplot = myfig5.add_subplot(111)
        aplot.plot(fit_x,fit_y, label="Data", linewidth=2)
        aplot.plot(fit_x, result.init_fit, 'k--', label="Initial Fit")
        aplot.plot(fit_x, result.best_fit, 'r--', label="Final Fit", linewidth=2)
        aplot.set_title("First Fit")
        print("Plotting components")
        for key in components:
            print(key)
            print(components[key])
            if key != "null_model" : 
                aplot.plot(fit_x, components[key], ':', label=key)
                aplot.legend(loc='best')
        #aplot.plot(fit_x, components[])
        canvas = FigureCanvasTkAgg(myfig5, f5)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, f5)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        return
    
    ##################################################################################################################################
    ##################################################################################################################################
    def baseline_all_fit(self):
        result, components, fit_x, fit_y = self.loaded_spectra.spectra[0].baseline_correct(self.reference_spectra, fit_bg=True, fit_scatter=True, fit_gaussian=False, wave_low=300, wave_high=800)
      
        resultwin = tk.Toplevel(width=600, height=800) #make new popup window
        resultwin.wm_title("Result of first fit")
      
        textbuffer = tk.Text(resultwin)
        textbuffer.insert(tk.END, result.fit_report())
        textbuffer.pack(side = tk.TOP, fill = tk.X, expand=1)
        f5 = tk.Frame(resultwin)
        f5.pack(side = tk.BOTTOM, fill = tk.BOTH, expand=1)

        myfig5 = plt.Figure(figsize=(5,5), dpi=100)
        aplot = myfig5.add_subplot(111)
        aplot.plot(fit_x,fit_y, label="Data", linewidth=2)
        aplot.plot(fit_x, result.init_fit, 'k--', label="Initial Fit")
        aplot.plot(fit_x, result.best_fit, 'r--', label="Final Fit", linewidth=2)
        aplot.set_title("First Fit")
        print("Plotting components")
        for key in components:
            print(key)
            print(components[key])
            if key != "null_model" : 
                aplot.plot(fit_x, components[key], ':', label=key)
                aplot.legend(loc='best')
        #aplot.plot(fit_x, components[])
        canvas = FigureCanvasTkAgg(myfig5, f5)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, f5)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        return




##################################################################################################################################
##################################################################################################################################
    def fitsimple(self):
      fit_titration_window = titrationwindow.TitrationWindow(master=self, spectra=self.loaded_spectra, reference_spectra=self.reference_spectra)
      return
##################################################################################################################################
##################################################################################################################################
    def fitcoop(self):
        return
##################################################################################################################################
##################################################################################################################################
    def about(self):
        tk.messagebox.showinfo("About this program", "P450 Titration Analysis by Donovan C. Haines, Sam Houston State University. Email:haines@shsu.edu.")
        return


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
class FitAllWindow(tk.Toplevel):
    def __init__(self,master=None, spectra=None, reference_spectra=None, width=600, height=800, fit_control_dict=None):
        tk.Toplevel.__init__(self, width=width, height=height)
        
        if fit_control_dict==None:
            self.fit_control_dict = { #default values
                "fit_bg": True,
                "fit_gaussian": False,
                "fit_scatter": False,
                "fit_ref1": True,
                "fit_ref2": True,
                "fit_ref3": False,
                "bg_val": 0.0,
                "scatter_int": 0.0,
                "scatter_int_fixed": False,
                "scatter_a": 1.0e4,
                "scatter_a-fixed": True,
                "scatter_power": 4.0,
                "scatter_power_fixed": False,
                "ref_coef1": 0.34,
                "ref_coef1_fixed": False,
                "ref_coef1_min": 0.0,
                "ref_coef2": 0.33,
                "ref_coef2_fixed": False,
                "ref_coef2_min": 0.0,
                "ref_coef3": 0.33,
                "ref_coef3_fixed": False,
                "ref_coef3_min": 0.0,
                "gaussian_int": 0.0,
                "gaussian_int_fixed": False,
                "gaussian_int_min": 0.0,
                "gaussian_cen": 260.0,
                "gaussian_cen_fixed": False,
                "gaussian_cen_min": 0.0,
                "guassian_cen_max": 1500.0,
                "gaussian_wid": 100.0,
                "guassian_wid_fixed": False,
                "gaussian-wid_min": 0.0,
                "guassian_wid_max": 1000.0,
                "wavelength_low": 0.0,
                "wavelength_high": 10000.0,
                "run_all": False}
        else:
            self.fit_control_dict=fit_control_dict
        
        self.master=master
        self.loaded_spectra=spectra
        self.reference_spectra=reference_spectra

        self.rowconfigure(0, weight=0, minsize=300)
        self.rowconfigure(1, weight=1, minsize=300)
        self.columnconfigure(0, weight=0, minsize=300)
        self.columnconfigure(1, weight=1, minsize=300)
        self.bv_use_bg = tk.BooleanVar()
        self.bv_use_bg.set(self.fit_control_dict.get("fit_bg", True))
        self.bv_use_gauss = tk.BooleanVar()
        self.bv_use_gauss.set(self.fit_control_dict.get("fit_gaussian", False))
        self.bv_use_scatter = tk.BooleanVar()
        self.bv_use_scatter.set(self.fit_control_dict.get("fit_scatter", False))
        self.bv_use_ref1 = tk.BooleanVar()
        self.bv_use_ref1.set(self.fit_control_dict.get("fit_ref1", True))
        self.bv_use_ref2 = tk.BooleanVar()
        self.bv_use_ref2.set(self.fit_control_dict.get("fit_ref2", True))
        self.bv_use_ref3 = tk.BooleanVar()
        self.bv_use_ref3.set(self.fit_control_dict.get("fit_ref3", True))
        self.select_var = tk.StringVar()
        self.bv_run_all = tk.BooleanVar()
        self.selected_spectrum = 0

        self.selection_frame = tk.Frame(self)
        self.selection_frame.grid(padx=30, pady=30)
     
        self.spectra_selector_opt=spectra.sample_names
        self.select_var.set(self.spectra_selector_opt[0])
        self.spectra_selector = tk.OptionMenu(self.selection_frame,self.select_var, *self.spectra_selector_opt)
        
        num_selection_columns=4
        self.spectra_selector.grid(columnspan=num_selection_columns)

        tk.Label(self.selection_frame, text="Constant?").grid(columnspan=2, row = 1, column = 3)

        self.b1 = tk.Checkbutton(self.selection_frame, text="Constant Background", variable = self.bv_use_bg)

        self.b1.grid(row=3, sticky="W")
        
        self.l1 = tk.Label(self.selection_frame, text="Initial constant bg value:")
        self.l1.grid(row=4, column=1, sticky="W")
        self.ev_bg_val = tk.DoubleVar()
        self.ent_bg_val = tk.Entry(self.selection_frame, textvariable=self.ev_bg_val)
        
        self.ent_bg_val.grid(row = 4, column = 2, sticky="W")
        self.bv_bg_val_fixed = tk.BooleanVar()
        self.b1a = tk.Checkbutton(self.selection_frame, variable = self.bv_bg_val_fixed)
        self.b1a.grid(row=4, column = 3, sticky="W", padx=20)


        self.b2 = tk.Checkbutton(self.selection_frame, text="Nearby Gaussian Background", variable = self.bv_use_gauss)
        self.b2.grid(row = 6, sticky = "W")
        
        self.l2a = tk.Label(self.selection_frame, text="Initial amplitude:")
        self.l2a.grid(row = 7, column=1, sticky = "W")
        self.ev_gaussian_amp = tk.DoubleVar()
        self.e2a = tk.Entry(self.selection_frame, textvariable = self.ev_gaussian_amp)
        self.e2a.grid(row = 7, column = 2, sticky = "W")
        self.bv_gaussian_amp_fixed = tk.BooleanVar()
        self.b2a = tk.Checkbutton(self.selection_frame, variable = self.bv_gaussian_amp_fixed)
        self.b2a.grid(row=7, column = 3, sticky="W", padx=20)

        self.l2b = tk.Label(self.selection_frame, text="Initial width:")
        self.l2b.grid(row = 8, column=1, sticky = "W")
        self.ev_gaussian_wid = tk.DoubleVar()
        self.e2b = tk.Entry(self.selection_frame, textvariable = self.ev_gaussian_wid)
        self.e2b.grid(row = 8, column = 2, sticky = "W")
        self.bv_gaussian_wid_fixed = tk.BooleanVar()
        self.b2b = tk.Checkbutton(self.selection_frame, variable = self.bv_gaussian_wid_fixed)
        self.b2b.grid(row=8, column = 3, sticky="W", padx=20)

        self.l2c = tk.Label(self.selection_frame, text="Initial wavelength:")
        self.l2c.grid(row = 9, column = 1, sticky = "W")
        self.ev_gaussian_cen = tk.DoubleVar()
        self.e2c = tk.Entry(self.selection_frame, textvariable = self.ev_gaussian_cen)
        self.e2c.grid(row = 9, column = 2, sticky = "W")
        self.bv_gaussian_cen_fixed = tk.BooleanVar()
        self.b2c = tk.Checkbutton(self.selection_frame, variable = self.bv_gaussian_cen_fixed)
        self.b2c.grid(row=9, column = 3, sticky="W", padx=20)

        self.b3 = tk.Checkbutton(self.selection_frame, text="Scattering Background", variable = self.bv_use_scatter)
        self.b3.grid(row=10, sticky="W")

        self.l3a = tk.Label(self.selection_frame, text="Initial amplitude:")
        self.l3a.grid(row = 11, column=1, sticky = "W")
        self.ev_scatter_int = tk.DoubleVar()
        self.e3a = tk.Entry(self.selection_frame, textvariable = self.ev_scatter_int)
        self.e3a.grid(row = 11, column = 2, sticky = "W")
        self.bv_scatter_int_fixed = tk.BooleanVar()
        self.b3a = tk.Checkbutton(self.selection_frame, variable = self.bv_scatter_int_fixed)
        self.b3a.grid(row=11, column = 3, sticky="W", padx=20)

        self.l3b = tk.Label(self.selection_frame, text="Initial proportionality:")
        self.l3b.grid(row = 12, column = 1, sticky = "W")
        self.ev_scatter_a = tk.DoubleVar()
        self.e3b = tk.Entry(self.selection_frame, textvariable = self.ev_scatter_a)
        self.e3b.grid(row = 12, column = 2, sticky = "W")
        self.bv_scatter_a_fixed = tk.BooleanVar()
        self.b3b = tk.Checkbutton(self.selection_frame, variable = self.bv_scatter_a_fixed)
        self.b3b.grid(row=12, column = 3, sticky="W", padx=20)

        self.l3c = tk.Label(self.selection_frame, text="Initial power:")
        self.l3c.grid(row = 13, column = 1, sticky = "W")
        self.ev_scatter_power = tk.DoubleVar()
        self.e3c = tk.Entry(self.selection_frame, textvariable = self.ev_scatter_power)
        self.e3c.grid(row = 13, column = 2, sticky = "W")
        self.bv_scatter_power_fixed = tk.BooleanVar()
        self.b3c = tk.Checkbutton(self.selection_frame, variable = self.bv_scatter_power_fixed)
        self.b3c.grid(row=13, column = 3, sticky="W", padx=20)

        self.b4 = tk.Checkbutton(self.selection_frame, text="Reference #1", variable = self.bv_use_ref1)
        self.b4.grid(row=14, sticky="W")

        self.l4a = tk.Label(self.selection_frame, text="Initial multiplier:")
        self.l4a.grid(row = 15, column=1, sticky = "W")
        self.ev_ref_coef1 = tk.DoubleVar()
        self.e4a = tk.Entry(self.selection_frame, textvariable = self.ev_ref_coef1)
        self.e4a.grid(row = 15, column = 2, sticky = "W")
        self.bv_ref_coef1_fixed = tk.BooleanVar()
        self.b4a = tk.Checkbutton(self.selection_frame, variable = self.bv_ref_coef1_fixed)
        self.b4a.grid(row=15, column = 3, sticky="W", padx=20)


        self.b5 = tk.Checkbutton(self.selection_frame, text="Reference #2", variable = self.bv_use_ref2)
        self.b5.grid(row=16, sticky="W")

        self.l5a = tk.Label(self.selection_frame, text="Initial multiplier:")
        self.l5a.grid(row = 17, column = 1, sticky = "W")
        self.ev_ref_coef2 = tk.DoubleVar()
        self.e5a = tk.Entry(self.selection_frame, textvariable = self.ev_ref_coef2)
        self.e5a.grid(row = 17, column = 2, sticky = "W")
        self.bv_ref_coef2_fixed = tk.BooleanVar()
        self.b5a = tk.Checkbutton(self.selection_frame, variable = self.bv_ref_coef2_fixed)
        self.b5a.grid(row=17, column = 3, sticky="W", padx=20)


        self.b6 = tk.Checkbutton(self.selection_frame, text="Reference #3", variable = self.bv_use_ref3)
        self.b6.grid(row = 18, sticky="W")
            
        self.l6a = tk.Label(self.selection_frame, text="Initial multiplier:")
        self.l6a.grid(row = 19, column=1, sticky = "W")
        self.ev_ref_coef3 = tk.DoubleVar()
        self.e6a = tk.Entry(self.selection_frame, textvariable = self.ev_ref_coef3)
        self.e6a.grid(row = 19, column = 2, sticky = "W")
        self.bv_ref_coef3_fixed = tk.BooleanVar()
        self.b7a = tk.Checkbutton(self.selection_frame, variable = self.bv_ref_coef3_fixed)
        self.b7a.grid(row=19, column = 3, sticky="W", padx=20)

        self.l8 = tk.Label(self.selection_frame, text="Fit between wavelengths ")
        self.l8.grid(row = 20, column = 0)
        self.ev_wavelength_low = tk.DoubleVar()
        self.ent_wave1 = tk.Entry(self.selection_frame, textvariable = self.ev_wavelength_low)
        self.ent_wave1.grid(row = 20, column = 1)
        self.l8a = tk.Label(self.selection_frame, text = " nm and ")
        self.l8a.grid(row = 20, column = 2)
        self.ev_wavelength_high = tk.DoubleVar()
        self.ent_wave2 = tk.Entry(self.selection_frame, textvariable = self.ev_wavelength_high)
        self.ent_wave2.grid(row = 20, column = 3)
        self.l8b = tk.Label(self.selection_frame, text = " nm")
        self.l8b.grid(row = 20, column = 4)
        self.ev_wavelength_low.set(self.fit_control_dict.get('wavelength_low', 0.0))
        self.ev_wavelength_high.set(self.fit_control_dict.get('wavelength_high', 1500.0))

        self.b_run_all = tk.Checkbutton(self.selection_frame, text="Use same parameters and fit all", variable = self.bv_run_all)
        self.b_run_all.grid(columnspan=num_selection_columns, row=24, padx=20, pady=30)
        
        self.b_last_fit = tk.Button(self.selection_frame, text="Set values to last fit result", command = self.set_to_last_fit)
        self.b_last_fit.grid(columnspan=num_selection_columns, row=25)

        self.b_ok = tk.Button(self.selection_frame, text="Run", command = self.baseline_fit_run)
        self.b_ok.grid(columnspan=num_selection_columns, row=26)

        self.selection_frame.grid(column = 0, row =0, padx=30, pady=30)
        self.update_fields_from_control(fit_control_dict)

        #########################################
        #Now define results and graph pane
        self.output_frame = tk.Frame(self, height=100)
        self.output_frame.grid(column = 0, row = 1, padx=20, pady=20, sticky='NSWE')
        self.results_scrollbar = tk.Scrollbar(self.output_frame)
        self.results_scrollbar.grid(column = 1, row = 0, sticky='NS')
        self.textbuffer = tk.Text(self.output_frame, yscrollcommand = self.results_scrollbar.set)
        self.textbuffer.grid(column = 0, row = 0, sticky='NSWE')
        self.results_scrollbar.config(command=self.textbuffer.yview)
        quote = "Results:"
        self.textbuffer.insert(tk.END, quote)

        self.graph_frame = tk.Frame(self, bg = "red", width=600, height=600)
        self.graph_frame.grid(rowspan=2, column = 1, row = 0, padx=20, pady=20, sticky='NSWE')

    ##############################################################################################################################
    ##############################################################################################################################
    def baseline_fit_run(self):
        if not self.bv_run_all.get():
            self.baseline_fit_run_each(self.get_spectrum_no())
            return
        for spec_no in range(0,self.loaded_spectra.number_spectra):
            self.baseline_fit_run_each(spec_no)
        return
        #TODO 3-31-19

    ##############################################################################################################################
    ##############################################################################################################################
    def baseline_fit_run_each(self, spectrum_no):
        this_control_dict = self.get_control_dict()
        print(this_control_dict)
        result, components, fit_x, fit_y, result_control = self.loaded_spectra.spectra[spectrum_no].baseline_correct(self.reference_spectra, control = this_control_dict)
        self.textbuffer.see("end")
        self.textbuffer.insert(tk.END, "\n")
        self.textbuffer.insert(tk.END, "*******************************************************\n")
        self.textbuffer.insert(tk.END, self.select_var.get())
        self.textbuffer.insert(tk.END, "\n*******************************************************\n")
        self.textbuffer.insert(tk.END, "\n")
        self.textbuffer.see("end")
        self.textbuffer.insert(tk.END, result.fit_report())
        self.textbuffer.see("end")
        #now update the graph and paramater windows, #TODO

        myfig5 = plt.Figure(figsize=(5,5), dpi=100)
        aplot = myfig5.add_subplot(111)
        aplot.plot(fit_x,fit_y, label="Data", linewidth=2)
        aplot.plot(fit_x, result.init_fit, 'k--', label="Initial Fit")
        aplot.plot(fit_x, result.best_fit, 'r--', label="Final Fit", linewidth=2)
        aplot.set_title("Fit of {0}".format(self.select_var.get()))
        print("Plotting components")
        for key in components:
            print(key)
            print(components[key])
            if key != "null_model" : 
                aplot.plot(fit_x, components[key], ':', label=key)
                aplot.legend(loc='best')
        
        #clear the self.graph_frame
        for widget in self.graph_frame.winfo_children():
            widget.destroy()

        #now set up new graph
        canvas = FigureCanvasTkAgg(myfig5, self.graph_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self.graph_frame)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


        return

    ##############################################################################################################################
    ##############################################################################################################################
    def set_to_last_fit(self):
        # pull the starting values from the results of the last fit if there is one
        if not hasattr(self.loaded_spectra.spectra[self.selected_spectrum], 'last_fit_result'):
            tk.messagebox.showerror("Error","I'm sorry. I was asked to retrieve the results of the last fit, but there does not appear to be a last fit to draw info from.")
            return

    ##############################################################################################################################
    ##############################################################################################################################
    def get_spectrum_no(self):
            if self.bv_run_all:
                return self.spectra_selector_opt.index(self.select_var.get()) #return list position of current selection
            else:
                return -1 #fit all so spectrum_no is returned as -1

    ##############################################################################################################################
    ##############################################################################################################################
    def get_control_dict(self):
        # return a control dictionary for the fit
        return_dict = { #default values
                "fit_bg": True,
                "fit_gaussian": False,
                "fit_scatter": False,
                "fit_ref1": True,
                "fit_ref2": True,
                "fit_ref3": True,
                "bg_val": 0.0,
                "scatter_int": 0.0,
                "scatter_int_fixed": False,
                "scatter_a": 1.0e4,
                "scatter_a-fixed": True,
                "scatter_power": 4.0,
                "scatter_power_fixed": False,
                "ref_coef1": 0.34,
                "ref_coef1_fixed": False,
                "ref_coef1_min": 0.0,
                "ref_coef2": 0.33,
                "ref_coef2_fixed": False,
                "ref_coef2_min": 0.0,
                "ref_coef3": 0.33,
                "ref_coef3_fixed": False,
                "ref_coef3_min": 0.0,
                "gaussian_int": 0.0,
                "gaussian_int_fixed": False,
                "gaussian_int_min": 0.0,
                "gaussian_cen": 260.0,
                "gaussian_cen_fixed": False,
                "gaussian_cen_min": 0.0,
                "guassian_cen_max": 1500.0,
                "gaussian_wid": 100.0,
                "guassian_wid_fixed": False,
                "gaussian-wid_min": 0.0,
                "guassian_wid_max": 1000.0,
                "wavelength_low": 0.0,
                "wavelength_high": 10000.0,
                "run_all": False}

        return_dict["fit_bg"] = self.bv_use_bg.get()
        return_dict["fit_gaussian"] = self.bv_use_gauss.get() 
        return_dict["fit_scatter"] = self.bv_use_scatter .get()
        return_dict["fit_ref1"] = self.bv_use_ref1.get()
        return_dict["fit_re2"] = self.bv_use_ref2.get()
        return_dict["fit_ref3"] = self.bv_use_ref3.get()
        #self.select_var 
        return_dict["run_all"] = self.bv_run_all.get()
        #constant background
        return_dict["bg_val"] = self.ev_bg_val.get()
        return_dict["bg_val_fixed"] = self.bv_bg_val_fixed.get()
        #scatter background
        return_dict["scatter_int"] = self.ev_scatter_int.get()
        return_dict["scatter_int_fixed"] = self.bv_scatter_int_fixed.get()
        return_dict["scatter_a"] = self.ev_scatter_a.get()
        return_dict["scatter_a_fixed"] = self.bv_scatter_a_fixed.get()
        return_dict["scatter_power"] = self.ev_scatter_power.get()
        return_dict["scatter_power_fixed"] = self.bv_scatter_power_fixed.get()
        #Gaussian Background
        return_dict["guassian_amp"] = self.ev_gaussian_amp.get()
        return_dict["guassian_amp_fixed"] = self.bv_gaussian_amp_fixed.get()
        return_dict["guassian_wid"] = self.ev_gaussian_wid.get()
        return_dict["guassian_wid_fixed"] = self.bv_gaussian_wid_fixed.get()
        return_dict["guassian_cen"] = self.ev_gaussian_cen.get()
        return_dict["guassian_cen_fixed"] = self.bv_gaussian_cen_fixed.get()
        return_dict["wavelength_low"] = self.ev_wavelength_low.get()
        return_dict["wavelength_high"] = self.ev_wavelength_high.get()
        return return_dict

    ##############################################################################################################################
    ##############################################################################################################################
    def update_fields_from_control(self, fit_con=None): 
        # Update the fit fields with the values in a fit control dictionary
        if fit_con == None:
            fit_con = {
                "empty_dict": True #to avoid None type erorrs when .get-ting from it.
                }
        
        default_dict = { #default values
                "fit_bg": True,
                "fit_gaussian": False,
                "fit_scatter": False,
                "fit_ref1": True,
                "fit_ref2": True,
                "fit_ref3": True,
                "bg_val": 0.0,
                "scatter_int": 0.0,
                "scatter_int_fixed": False,
                "scatter_a": 1.0e4,
                "scatter_a-fixed": True,
                "scatter_power": 4.0,
                "scatter_power_fixed": False,
                "ref_coef1": 0.34,
                "ref_coef1_fixed": False,
                "ref_coef1_min": 0.0,
                "ref_coef2": 0.33,
                "ref_coef2_fixed": False,
                "ref_coef2_min": 0.0,
                "ref_coef3": 0.33,
                "ref_coef3_fixed": False,
                "ref_coef3_min": 0.0,
                "gaussian_int": 0.0,
                "gaussian_int_fixed": False,
                "gaussian_int_min": 0.0,
                "gaussian_cen": 260.0,
                "gaussian_cen_fixed": False,
                "gaussian_cen_min": 0.0,
                "guassian_cen_max": 1500.0,
                "gaussian_wid": 100.0,
                "guassian_wid_fixed": False,
                "gaussian-wid_min": 0.0,
                "guassian_wid_max": 1000.0,
                "wavelength_low": 0.0,
                "wavelength_high": 10000.0,
                "run_all": False}
        self.bv_use_bg.set(fit_con.get("fit_bg", default_dict.get("fit_bg", True)))
        self.bv_use_gauss.set(fit_con.get("fit_gaussian", default_dict.get("fit_gaussian", False))) 
        self.bv_use_scatter.set(fit_con.get("fit_scatter", default_dict.get("fit_scatter", False)))
        self.bv_use_ref1.set(fit_con.get("fit_ref1", default_dict.get("fit_ref1", True)))
        self.bv_use_ref2.set(fit_con.get("fit_ref2", default_dict.get("fit_ref2", True)))
        self.bv_use_ref3.set(fit_con.get("fit_ref3", default_dict.get("fit_ref3", False)))
        #self.select_var 
        self.bv_run_all.set(fit_con.get("run_all", default_dict.get("run_all", True)))
        #constant background
        self.ev_bg_val.set(fit_con.get("bg_val", default_dict.get("bg_val", 0.0)))
        self.bv_bg_val_fixed.set(fit_con.get("bg_val_fixed", default_dict.get("bg_val_fixed", False)))
        #scatter background
        self.ev_scatter_int.set(fit_con.get("scatter_int", default_dict.get("scatter_int", 0.0)))
        self.bv_scatter_int_fixed.set(fit_con.get("scatter_int_fixed", default_dict.get("scatter_int_fixed", 0.0)))
        self.ev_scatter_a.set(fit_con.get("scatter_a", default_dict.get("scatter_a", 1e4)))
        self.bv_scatter_a_fixed.set(fit_con.get("scatter_a_fixed", default_dict.get("scatter_a_fixed", True)))
        self.ev_scatter_power.set(fit_con.get("scatter_power", default_dict.get("scatter_power", 4.0)))
        self.bv_scatter_power_fixed.set(fit_con.get("scatter_power_fixed", default_dict.get("scatter_power_fixed", False)))
        #Gaussian Background
        self.ev_gaussian_amp.set(fit_con.get("guassian_amp", default_dict.get("guassian_amp", 0.0)))
        self.bv_gaussian_amp_fixed.set(fit_con.get("guassian_amp_fixed", default_dict.get("guassian_amp_fixed", False)))
        self.ev_gaussian_wid.set(fit_con.get("guassian_wid", default_dict.get("guassian_wid", 50.0)))
        self.bv_gaussian_wid_fixed.set(fit_con.get("guassian_wid_fixed", default_dict.get("guassian_wid_fixed", False)))
        self.ev_gaussian_cen.set(fit_con.get("guassian_cen", default_dict.get("guassian_cen", 0.0)))
        self.bv_gaussian_cen_fixed.set(fit_con.get("guassian_cen_fixed", default_dict.get("guassian_cen_fixed", False)))
        self.ev_wavelength_low.set(fit_con.get("wavelength_low", default_dict.get("wavelength_low", 0.0)))
        self.ev_wavelength_high.set(fit_con.get("wavelength_high", default_dict.get("wavelength_high", 0.0)))
        return 

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
# Main

if __name__ == '__main__':
    #root = tk.Tk()
    my_gui = MyGUI()
    
    my_gui.mainloop()



