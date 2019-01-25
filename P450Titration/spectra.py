import pandas
from tkinter import *

import numpy as np
from scipy.optimize import optimize

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
class uvvis_spectrum:
  def __init__(self):
    self.wl_low=-1.0
    self.wl_high=-1.0
    self.filename=""
    self.sample_name=""
    self.collection_time=""
    self.instrument=""
    self.avg_time=0.0
    self.baseline_method=""
    self.comments=""
    self.wavelength=list()
    self.absorbance=list()
  
    ####################################################################
  # borrowed from https://stackoverflow.com/questions/1535327/how-to-print-objects-of-class-using-print
  def __str__(self):
    return str(self.__class__) + ": " + str(self.__dict__)   
  

def baseline_correct(self, reference_spectra=None, fit_bg=True, fit_scatter=True, fit_gaussian=True):
    if reference_spectra==None:
        numref=0
        def fit_function(wavelength, )
    else:
        numref = reference_spectra.number_spectra
    


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
class spectral_collection:
  def __init__(self, filename="", load=FALSE):
    self.number_spectra=0
    self.number_files_loaded=0
    self.spectra=list()
    self.filenames=list()
    self.rawdata=list()
    self.df = pandas.DataFrame()
    if load : self.load(filename)
  
  ####################################################################
  # borrowed from https://stackoverflow.com/questions/1535327/how-to-print-objects-of-class-using-print
  def __str__(self):
    return str(self.__class__) + ": " + str(self.__dict__) 
    
  ####################################################################
  def load(self, filename="", generatedf=True):
    #load a csv file from the Cary50
    if filename=="": 
      filename = filedialog.askopenfilename()
    print(filename)
    
    #Add filename to list of filenames (because can import more files later)
    self.filenames.append(filename)
    self.rawdata.append(pandas.read_csv(filename, header=None))
    self.number_files_loaded = self.number_files_loaded + 1


    names=self.rawdata[self.number_files_loaded-1].iloc[0]
    labels=self.rawdata[self.number_files_loaded-1].iloc[1]
    print(names)
    print(labels)

    # Determine number of valid wavelength, Abs rows
    num_wavelength=0
    for count in range(self.rawdata[self.number_files_loaded-1].shape[0]):
      try: #test first to columns to see if they are floats
        #print("Cycle {}".format(count))
        #print(self.rawdata[self.number_files_loaded-1].loc[count,0:1])
        float(self.rawdata[self.number_files_loaded-1].loc[count,0])
        float(self.rawdata[self.number_files_loaded-1].loc[count,1])
        num_wavelength = num_wavelength + 1
        #print("Count {} : yes".format(count))
      except:
        print("Count {} : not a float pair row.".format(count))
    print("Found {} rows of float pair data.".format(num_wavelength))

    # for each spectrum
    for count in range((self.rawdata[self.number_files_loaded-1].shape[1])//2):
      #print("Spectral loop count = {}".format(count))
      this_spectrum = uvvis_spectrum()
      this_spectrum.filename = filename
      #this_spectrum.wl_low = min()
      xvalues=pandas.to_numeric(self.rawdata[self.number_files_loaded-1].loc[2:num_wavelength+1, count*2])
      yvalues=pandas.to_numeric(self.rawdata[self.number_files_loaded-1].loc[2:num_wavelength+1, (count*2)+1])
      this_spectrum.wavelengths = xvalues
      this_spectrum.absorbance = yvalues
      this_spectrum.sample_name=(self.rawdata[self.number_files_loaded-1].loc[0, (count*2)+1])

      self.spectra.append(this_spectrum)
      self.number_spectra = self.number_spectra+1

      if generatedf:
        xvalues.name="WL:"+names[count*2]
        yvalues.name="AU:"+names[count*2]
        #print(type(xvalues))
        
        tempdf = pandas.concat([xvalues,yvalues], axis=1, ignore_index=True)
        tempdf.columns=[xvalues.name, yvalues.name] # not sure why this is necessary, but seems to be
        #print(tempdf)
        if count==0:
          self.df=tempdf
        else:
          self.df = pandas.concat([self.df,tempdf], axis=1)
          #print("Appended df")
          #print("xvalues")
          #print(xvalues[0:3])
          #print("yvalues")
          #print(yvalues[0:3])
          #print("tempdf")
          #print(tempdf.iloc[0:3,:])
          #print("self.df:")
          #print(self.df.iloc[0:3,:])
        print("self.df in spectral collection ended up being:")
        print(self.df)
