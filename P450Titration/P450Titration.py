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

class MyGUI:
    def __init__(self,master):
      self.master = master
      master.title("P450 Titration Analysis")
      master.geometry('1920x1080')
      self.loaded_spectra = None
      self.reference_spectra = None
      
      #self.greet_button = Button(master, text="Greet", command=self.greet)
      #self.greet_button.pack()
      
      #self.close_button = Button(master, text="Close", command=master.destroy)
      #self.close_button.pack()

      menu=tk.Menu(master)
      new_item=tk.Menu(menu)
      new_item.add_command(label="Load Spectra", command=self.load_spectra)
      new_item.add_command(label="Load Reference Spectra", command=self.load_reference)
      new_item.add_command(label="Exit", command=master.destroy)
      menu.add_cascade(label='File', menu=new_item)
      master.config(menu=menu)
      
      
      new_item2=tk.Menu(menu)
      new_item2.add_command(label="Baseline Correct", command=self.baseline)
      new_item2.add_command(label="Baseline Correct All", command=self.baseline_all)
      new_item2.add_command(label="Fit simple binding", command=self.fitsimple)
      new_item2.add_command(label="Fit coop binding", command=self.fitcoop)
      menu.add_cascade(label="Analyze", menu=new_item2)

      new_item3=tk.Menu(menu)
      new_item3.add_command(label="About this program", command=self.about)
      menu.add_cascade(label="Help", menu=new_item3)

          
    ################################################################################################
    def load_spectra(self):
      self.loaded_spectra = loaded_spectra = spectra.spectral_collection(load=True)
      print(loaded_spectra)

      f1=tk.Frame(self.master, width=400)
      f1.pack(side = tk.LEFT, fill=tk.Y,expand=0)
      df = loaded_spectra.df
      pt = pdt.Table(f1, dataframe=df,
                                    showtoolbar=True, showstatusbar=True)
      pt.show()

      f2 = tk.Frame(self.master)
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


    def load_reference(self):
      #self.ref_filename = filedialog.askopenfilename()
      #print(self.ref_filename)

      self.reference_spectra = reference_spectra = spectra.spectral_collection(load=True)
      #self.ref_data=pd.read_csv(self.ref_filename, header=None)
      #f1=Frame(self.master)
      #f1.pack(side = LEFT, fill=BOTH,expand=1)
      df = reference_spectra.df
      #pt = Table(f1, dataframe=df,
      #                              showtoolbar=True, showstatusbar=True)
      #pt.show()

      f2 = tk.Frame(self.master)
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
    
    ############################################################################
    def baseline_all(self):
        base_win=tk.Toplevel(width=600, height=800) #make new popup window
        base_win.rowconfigure(0, weight=0, minsize=300)
        base_win.rowconfigure(1, weight=1, minsize=300)
        base_win.columnconfigure(0, weight=0, minsize=300)
        base_win.columnconfigure(1, weight=1, minsize=300)
        bv_use_bg = tk.IntVar()
        bv_use_gauss = tk.IntVar()
        bv_use_scatter = tk.IntVar()
        bv_use_ref1 = tk.IntVar()
        bv_use_ref2 = tk.IntVar()
        bv_use_ref3 = tk.IntVar()
        select_var = tk.StringVar()
        bv_run_all = tk.IntVar()

        selection_frame = tk.Frame(base_win)
        selection_frame.grid(padx=30, pady=30)
        #sample_names=list()
        #for this_spectrum in self.loaded_spectra.spectra:
        #    sample_names.append(this_spectrum.sample_name)
        
        #Create dropdown box listing all spectra available (non-ref)
        w_opt=self.loaded_spectra.sample_names
        select_var.set(w_opt[0])
        w = tk.OptionMenu(selection_frame,select_var, *w_opt)
        #w = OptionMenu(selection_frame, select_var, self.loaded_spectra.sample_names)
        num_selection_columns=4
        w.grid(columnspan=num_selection_columns)

        tk.Label(selection_frame, text="Constant?").grid(columnspan=2, row = 1, column = 3)

        b1 = tk.Checkbutton(selection_frame, text="Constant Background", variable = bv_use_bg)
        b1.grid(row=3, sticky="W")
        
        l1 = tk.Label(selection_frame, text="Initial constant bg value:")
        l1.grid(row=4, column=1, sticky="W")
        e1 = tk.Entry(selection_frame)
        e1.grid(row = 4, column = 2, sticky="W")
        bv_b1a = tk.DoubleVar()
        b1a = tk.Checkbutton(selection_frame, variable = bv_b1a)
        b1a.grid(row=4, column = 3, sticky="W", padx=20)


        b2 = tk.Checkbutton(selection_frame, text="Nearby Gaussian Background", variable = bv_use_gauss)
        b2.grid(row = 6, sticky = "W")
        
        l2a = tk.Label(selection_frame, text="Initial amplitude:")
        l2a.grid(row = 7, column=1, sticky = "W")
        e2a = tk.Entry(selection_frame)
        e2a.grid(row = 7, column = 2, sticky = "W")
        bv_b2a = tk.DoubleVar()
        b2a = tk.Checkbutton(selection_frame, variable = bv_b2a)
        b2a.grid(row=7, column = 3, sticky="W", padx=20)

        l2b = tk.Label(selection_frame, text="Initial width:")
        l2b.grid(row = 8, column=1, sticky = "W")
        e2b = tk.Entry(selection_frame)
        e2b.grid(row = 8, column = 2, sticky = "W")
        bv_b2b = tk.DoubleVar()
        b2b = tk.Checkbutton(selection_frame, variable = bv_b2b)
        b2b.grid(row=8, column = 3, sticky="W", padx=20)

        l2c = tk.Label(selection_frame, text="Initial wavelength:")
        l2c.grid(row = 9, column = 1, sticky = "W")
        e2c = tk.Entry(selection_frame)
        e2c.grid(row = 9, column = 2, sticky = "W")
        bv_b2c = tk.DoubleVar()
        b2c = tk.Checkbutton(selection_frame, variable = bv_b2c)
        b2c.grid(row=9, column = 3, sticky="W", padx=20)

        b3 = tk.Checkbutton(selection_frame, text="Scattering Background", variable = bv_use_scatter)
        b3.grid(row=10, sticky="W")

        l3a = tk.Label(selection_frame, text="Initial amplitude:")
        l3a.grid(row = 11, column=1, sticky = "W")
        e3a = tk.Entry(selection_frame)
        e3a.grid(row = 11, column = 2, sticky = "W")
        bv_b3a = tk.DoubleVar()
        b3a = tk.Checkbutton(selection_frame, variable = bv_b3a)
        b3a.grid(row=11, column = 3, sticky="W", padx=20)


        l3b = tk.Label(selection_frame, text="Initial proportionality:")
        l3b.grid(row = 12, column = 1, sticky = "W")
        e3b = tk.Entry(selection_frame)
        e3b.grid(row = 12, column = 2, sticky = "W")
        bv_b3b = tk.DoubleVar()
        b3b = tk.Checkbutton(selection_frame, variable = bv_b3b)
        b3b.grid(row=12, column = 3, sticky="W", padx=20)

        l3c = tk.Label(selection_frame, text="Initial power:")
        l3c.grid(row = 13, column = 1, sticky = "W")
        e3c = tk.Entry(selection_frame)
        e3c.grid(row = 13, column = 2, sticky = "W")
        bv_b3c = tk.DoubleVar()
        b3c = tk.Checkbutton(selection_frame, variable = bv_b3c)
        b3c.grid(row=13, column = 3, sticky="W", padx=20)

        b4 = tk.Checkbutton(selection_frame, text="Reference #1", variable = bv_use_ref1)
        b4.grid(row=14, sticky="W")

        l4a = tk.Label(selection_frame, text="Initial multiplier:")
        l4a.grid(row = 15, column=1, sticky = "W")
        e4a = tk.Entry(selection_frame)
        e4a.grid(row = 15, column = 2, sticky = "W")
        bv_b4a = tk.DoubleVar()
        b4a = tk.Checkbutton(selection_frame, variable = bv_b4a)
        b4a.grid(row=15, column = 3, sticky="W", padx=20)


        b5 = tk.Checkbutton(selection_frame, text="Reference #2", variable = bv_use_ref2)
        b5.grid(row=16, sticky="W")

        l5a = tk.Label(selection_frame, text="Initial multiplier:")
        l5a.grid(row = 17, column = 1, sticky = "W")
        e5a = tk.Entry(selection_frame)
        e5a.grid(row = 17, column = 2, sticky = "W")
        bv_b5a = tk.DoubleVar()
        b5a = tk.Checkbutton(selection_frame, variable = bv_b5a)
        b5a.grid(row=17, column = 3, sticky="W", padx=20)


        b6 = tk.Checkbutton(selection_frame, text="Reference #3", variable = bv_use_ref3)
        b6.grid(row = 18, sticky="W")
            
        l6a = tk.Label(selection_frame, text="Initial multiplier:")
        l6a.grid(row = 19, column=1, sticky = "W")
        e6a = tk.Entry(selection_frame)
        e6a.grid(row = 19, column = 2, sticky = "W")
        bv_b7a = tk.DoubleVar()
        b7a = tk.Checkbutton(selection_frame, variable = bv_b7a)
        b7a.grid(row=19, column = 3, sticky="W", padx=20)

        l8 = tk.Label(selection_frame, text="Fit between wavelengths ")
        l8.grid(row = 20, column = 0)
        ent_wave1 = tk.Entry(selection_frame)
        ent_wave1.grid(row = 20, column = 1)
        l8a = tk.Label(selection_frame, text = " nm and ")
        l8a.grid(row = 20, column = 2)
        ent_wave2 = tk.Entry(selection_frame)
        ent_wave2.grid(row = 20, column = 3)
        l8b = tk.Label(selection_frame, text = " nm")
        l8b.grid(row = 20, column = 4)

        b_run_all = tk.Checkbutton(selection_frame, text="Use same parameters and fit all", variable = bv_run_all)
        b_run_all.grid(columnspan=num_selection_columns, row=25, padx=20, pady=30)

        b_ok = tk.Button(selection_frame, text="Run", command = lambda: self.baseline())
        b_ok.grid(columnspan=num_selection_columns, row=26)

        selection_frame.grid(column = 0, row =0, padx=30, pady=30)
        
        # now set to not shrink window too small for entry widgets
        #base_win.update()
        #temp=selection_frame.grid_info()
        #rowconfigure(0, minsize=selection_frame.winfo_height)
        #base_win.columnconfigure(0, minsize=selection_frame.winfo_width)


        #########################################
        #Now define results and graph pane
        output_frame = tk.Frame(base_win, height=100)
        output_frame.grid(column = 0, row = 1, padx=20, pady=20, sticky='NSWE')
        results_scrollbar = tk.Scrollbar(output_frame)
        results_scrollbar.grid(column = 1, row = 0, sticky='NS')
        textbuffer = tk.Text(output_frame, yscrollcommand = results_scrollbar.set)
        textbuffer.grid(column = 0, row = 0, sticky='NSWE')
        results_scrollbar.config(command=textbuffer.yview)
        quote = """HAMLET: To be, or not to be--that is the question:
        Whether 'tis nobler in the mind to suffer
        The slings and arrows of outrageous fortune
        Or to take arms against a sea of troubles
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .
        .

        And by opposing end them. To die, to sleep--
        No more--and by a sleep to say we end
        The heartache, and the thousand natural shocks
        That flesh is heir to. 'Tis a consummation
        Devoutly to be wished."""
        textbuffer.insert(tk.END, quote)

        graph_frame = tk.Frame(base_win, bg = "red", width=600, height=600)
        graph_frame.grid(rowspan=2, column = 1, row = 0, padx=20, pady=20, sticky='NSWE')


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

    def fitsimple(self):
   
      return

    def fitcoop(self):
        return

    def about(self):
        tk.messagebox.showinfo("About this program", "P450 Titration Analysis by Donovan C. Haines, Sam Houston State University. Email:haines@shsu.edu.")
        return



root = tk.Tk()
my_gui = MyGUI(root)
root.mainloop()



