from tkinter import *  #to use graphical user interface (windowed program)
from tkinter import filedialog
import tkinter.messagebox
import pandas
from pandastable import Table, TableModel

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

      
      #self.greet_button = Button(master, text="Greet", command=self.greet)
      #self.greet_button.pack()
      
      #self.close_button = Button(master, text="Close", command=master.destroy)
      #self.close_button.pack()

      menu=Menu(master)
      new_item=Menu(menu)
      new_item.add_command(label="Load Spectra", command=self.load_spectra)
      new_item.add_command(label="Load Reference Spectra", command=self.load_reference)
      new_item.add_command(label="Exit", command=master.destroy)
      menu.add_cascade(label='File', menu=new_item)
      master.config(menu=menu)
      
      
      new_item2=Menu(menu)
      new_item2.add_command(label="Baseline Correct", command=self.baseline)
      new_item2.add_command(label="Fit simple binding", command=self.fitsimple)
      new_item2.add_command(label="Fit coop binding", command=self.fitcoop)
      menu.add_cascade(label="Analyze", menu=new_item2)

      new_item3=Menu(menu)
      new_item3.add_command(label="About this program", command=self.about)
      menu.add_cascade(label="Help", menu=new_item3)


    
    def greet(self):
      print("Greetings!")
    
    def load_spectra(self):
      loaded_spectra = spectra.spectral_collection(load=TRUE)
      print(loaded_spectra)

      f1=Frame(self.master)
      f1.pack(side = LEFT, fill=BOTH,expand=1)
      df = loaded_spectra.df
      pt = Table(f1, dataframe=df,
                                    showtoolbar=True, showstatusbar=True)
      pt.show()

      f2 = Frame(self.master)
      f2.pack(side = RIGHT, fill = BOTH, expand=1)

      myfig = plt.Figure(figsize=(5,5), dpi=100)
      aplot = myfig.add_subplot(111)
      aplot.plot(loaded_spectra.spectra[0].wavelengths,loaded_spectra.spectra[0].absorbance)
      for count in range(1,loaded_spectra.number_spectra):
          aplot.plot(loaded_spectra.spectra[count].wavelengths,loaded_spectra.spectra[count].absorbance)

      canvas = FigureCanvasTkAgg(myfig, f2)
      canvas.draw()
      canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=TRUE)
      toolbar = NavigationToolbar2Tk(canvas, f2)
      toolbar.update()
      canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=TRUE)


    def load_reference(self):
      self.ref_filename = filedialog.askopenfilename()
      print(self.ref_filename)
      self.ref_data=pandas.read_csv(self.ref_filename, header=None)
    
    def baseline(self):
        return

    def fitsimple(self):
        return

    def fitcoop(self):
        return

    def about(self):
        tkinter.messagebox.showinfo("About this program", "P450 Titration Analysis by Donovan C. Haines, Sam Houston State University")
        return



root = Tk()
my_gui = MyGUI(root)
root.mainloop()



