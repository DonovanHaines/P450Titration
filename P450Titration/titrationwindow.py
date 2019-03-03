import tkinter as tk  #to use graphical user interface (windowed program)
from tkinter import filedialog
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
    def __init__(self,master=None, spectra=None, reference_spectra=None, width=None, height=None, fit_control_dict=None):
        
        if master==None:
           if width==None:
               width=800
           if height==None:
               height=800
        else:
            if width==None:
                width = int(master.winfo_width()*0.9)
            if height==None:
                height= int(master.winfo_height()*0.9)
        print("Titration window width {0} and height {1} from locals".format(width,height))
        
        tk.Toplevel.__init__(self, width=width, height=height)
        self.geometry="{0}x{1}".format(width, height)
        print("Titration window width {0} and height {1} from self.winfo_".format(self.winfo_width(),self.winfo_height()))
        self.master=master
        self.loaded_spectra=spectra
        self.reference_spectra=reference_spectra
    
        #   Set up the menu system
        menu=tk.Menu(self)
        new_item=tk.Menu(menu)
        new_item.add_command(label="Load Volume Data", command=self.load_volumes)
        new_item.add_command(label="Save Table", command=self.save_table)
        new_item.add_command(label="Exit", command=master.destroy)
        menu.add_cascade(label='File', menu=new_item)
        self.config(menu=menu)
     
        new_item2=tk.Menu(menu)
        new_item2.add_command(label="Simple Fit", command=self.fit_simple)
        new_item2.add_command(label="Cooperative Fit", command=self.fit_coop)

        menu.add_cascade(label="Fit", menu=new_item2)
     
        new_item3=tk.Menu(menu)
        new_item3.add_command(label="About this program", command=self.about)
        menu.add_cascade(label="Help", menu=new_item3)
    
        self.grid_columnconfigure(0,weight=1)
        self.grid_rowconfigure(2, weight=1)

        # Now to add controls
        self.control_frame = tk.Frame(self)
        self.control_frame.grid(row = 0, column = 0, padx=10, pady=5)
        self.l_fill=tk.Label(self.control_frame, text="Fill Column with Value:")
        self.fill_value = tk.StringVar()
        column_names = ['Name', 'Include','Total_Volume', 'Other_Volume', 'Enzyme_Volume', 'Substrate_Volume',
                        'Sub_Stock','Enz_Stock', 'S_Total', 'Enz_Total']
        data = pd.DataFrame(np.nan, index=spectra.sample_names, columns=column_names)
        self.ev_fill = tk.Entry(self.control_frame, textvariable=self.fill_value)
        self.selected_column = tk.StringVar()
        self.selected_column.set(column_names[0])
        self.column_selector = tk.OptionMenu(self.control_frame, self.selected_column, *column_names)
        self.bv_fill = tk.BooleanVar()
        self.b_fill = tk.Button(self.control_frame, text="Fill", command=self.column_fill)
        self.l_fill.grid(column=0, row=0)
        self.ev_fill.grid(column=1, row=0)
        self.column_selector.grid(column = 2, row = 0)
        self.b_fill.grid(column = 3, row = 0)

        #add columns
        add_column_choices = ['Absorbance', 'Absorbance Difference', 'Reference', 
                              'Ref Fraction', 'Baseline Parameter']
        self.l_addcol=tk.Label(self.control_frame, text="Add Column of type:")
        self.l_addcol.grid(column = 0, row = 1)
        self.selected_add_column_choice = tk.StringVar()
        self.selected_add_column_choice.set(add_column_choices[0])
        self.add_column_selector = tk.OptionMenu(self.control_frame, 
                                                 self.selected_add_column_choice, 
                                                 *add_column_choices)
        self.bv_add_col = tk.BooleanVar()
        self.b_add_col = tk.Button(self.control_frame, text = "Add Column", command=self.add_column)
        self.add_column_selector.grid(row=1, column=1, padx=(0,0))
        self.b_add_col.grid(row = 1, column = 2)

        self.text_frame = tk.Frame(self, background="white")
        self.l_click=tk.Label(self.text_frame, text="Double click cells to edit individually.")
        self.l_click.grid(column = 0, row = 1)
        self.text_frame.grid(row=1, column=0, pady=10)

        # Now to add the table
        print("Titration window width {0} and height {1} at spreadsheet frame formation".format(self.winfo_width(),self.winfo_height()))
        self.spreadsheet_frame = tk.Frame(self,  height=self.winfo_height(), width=self.winfo_width(), background="red")
        self.spreadsheet_frame.grid(column = 0, row = 2, sticky='nswe', columnspan=100)

        rows = spectra.number_spectra
        

        #self.data = pd.DataFrame(np.nan, index=spectra.sample_names, columns=column_names)
        print("titrationwindow data is currently:")
        print(data)
        self.pt = pdt.Table(self.spreadsheet_frame, dataframe=data, showtoolbar=True, showstatusbar=True)
        self.pt.model.df['Name']=list(self.pt.model.df.index)
        self.pt.model.df['Include']='True'
        self.pt.model.df['Total_Volume']="1000.0"
        self.pt.show()
        self.pt.redraw()
        print("data has been displayed.")
        self.geometry="{0}x{1}".format(width, height)
        print("Titration window width {0} and height {1} from self.winfo_".format(self.winfo_width(),self.winfo_height()))
    
    def add_column(self):
        print(self.selected_add_column_choice.get())
        switcher_choices = {'Absorbance': self.add_col_abs, 
                            'Absorbance Difference': self.add_col_diff, 
                            'Reference': self.add_col_ref, 
                            'Ref Fraction': self.add_col_reffrac, 
                            'Baseline Parameter': self.add_col_baseline}
        function_to_call = switcher_choices.get(self.selected_add_column_choice.get(), "Invalid column choice")
        print("Add column function to call is")
        print(function_to_call)
        function_to_call()
        return
    
    def add_col_abs(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def add_col_diff(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def add_col_ref(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def add_col_reffrac(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def add_col_baseline(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return


    def refresh_column_selector(self):
        #  Reset var and delete all old options
        self.selected_column.set('')
        self.column_selector['menu'].delete(0, 'end')

        # Insert list of new options (tk._setit hooks them up to self.selected_column)
        new_choices = self.pt.model.df.columns
        for choice in new_choices:
            self.column_selector['menu'].add_command(label=choice, command=tk._setit(self.selected_column, choice))
        
    def load_volumes(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def save_table(self):
        f = tk.filedialog.asksaveasfile(mode='w', defaultextension=".txt")
        if f is None: 
            return
        self.pt.model.df.to_csv(f)
        return

    def fit_simple(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def fit_coop(self):
        tk.messagebox.showinfo("OOps!", "My apologies, this is not yet implemented.")
        return

    def about(self):
        tk.messagebox.showinfo("About this program", "P450 Titration Analysis by Donovan C. Haines, Sam Houston State University. Email:haines@shsu.edu.")
        return

    def column_fill(self):
        self.pt.model.df[self.selected_column.get()] = self.fill_value.get()
        self.pt.redraw()
        print("column filled:")
        print(self.selected_column.get())
        print(self.pt.model.df)