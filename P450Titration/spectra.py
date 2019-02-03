import pandas as pd
import tkinter as tk
from tkinter import filedialog
import numpy as np
import math

from scipy.optimize import optimize
from lmfit import Model, report_fit, report_errors

import bisect

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
    self.wavelength=[]
    self.absorbance=[]
    self.last_fit_result=None
    
  
    ####################################################################
  # borrowed from https://stackoverflow.com/questions/1535327/how-to-print-objects-of-class-using-print
  def __str__(self):
    return str(self.__class__) + ": " + str(self.__dict__)   
  
  def last_fit_was_gaussian():
    if self.last_fit_result==None: 
        return False
    

  def baseline_correct(self, reference_spectra=None, fit_bg=True, fit_scatter=True, fit_gaussian=True, wave_low=300, wave_high=800, hold_scatter_power=True):
    if reference_spectra==None:
        numref=0
        print("Numref")
        print(numref)
    else:
      if type(reference_spectra==spectral_collection):
        numref = reference_spectra.number_spectra
      else:
        print("Warning - Exiting baseline_correct due to wrong reference_spectra type.")
        return

    mask = (wave_low <= self.wavelength) & (self.wavelength <= wave_high)
    fit_x = np.asarray(self.wavelength[mask])
    fit_y = np.asarray(self.absorbance[mask])


    if numref >= 1:
        ref1seq = np.empty(len(fit_x))
        for count in range(len(fit_x)):
            ref1seq[count] = np.abs(fit_x[count]-reference_spectra.spectra[0].wavelength[mask]).idxmin() #precalculate indices for each wavelength
    if numref >= 2:
        ref2seq = np.empty(len(fit_x))
        for count in range(len(fit_x)):
            ref2seq[count] = np.abs(fit_x[count]-reference_spectra.spectra[1].wavelength[mask]).idxmin() #precalculate indices for each wavelength
    if numref >= 3:
        ref3seq = np.empty(len(fit_x))
        for count in range(len(fit_x)):
            ref3seq[count] = np.abs(fit_x[count]-reference_spectra.spectra[2].wavelength[mask]).idxmin() #precalculate indices for each wavelength
 

    #####model functions
    def bg(x, bg_val):
        return np.asarray([bg_val] * len(x))
   
    def scatter(x, scatter_int, scatter_a, scatter_power):  #or fix scatter_power to -4
        return np.asarray(scatter_int*np.log10(1/(1-scatter_a*(x**scatter_power))))
    
    def ref1(x, refcoef1):
        return np.asarray(refcoef1*reference_spectra.spectra[0].absorbance[ref1seq])

    def ref2(x, refcoef2):
        return np.asarray(refcoef2*reference_spectra.spectra[1].absorbance[ref2seq])
    
    def ref3(x, refcoef3):
        return np.asarray(refcoef3*reference_spectra.spectra[1].absorbance[ref3seq])
    
    #def ref2(x, refcoef1, refcoef2):
    #   
    #    return refcoef1*reference_spectra.spectra[0].absorbance[ref1seq] + \
    #           refcoef2*reference_spectra.spectra[1].absorbance[ref2seq]
    #     
   # 
   # def ref3(x, refcoef1, refcoef2, refcoef3):
   #     return refcoef1*reference_spectra.spectra[0].absorbance[ref1seq] + \
   #            refcoef2*reference_spectra.spectra[1].absorbance[ref2seq] + \
   #            refcoef3*reference_spectra.spectra[2].absorbance[ref2seq]
    
    def gaussian(x, amp, cen, wid):
      #"1-d gaussian: gaussian(x, amp, cen, wid)"
      return np.asarray([amp * np.exp(-(xt-cen)**2.0 / wid) for xt in x])
     
    def null_model(x):
        return np.zeros(len(x))


    if fit_gaussian: 
        gaussian_part = Model(gaussian)
    else: 
        gaussian_part = Model(null_model)

    if fit_scatter:
        scatter_part = Model(scatter)
    else:
        scatter_part = Model(null_model)
    
    if fit_bg:
        bg_part = Model(bg)
    else:
        bg_part = Model(null_model)

    if numref>0: reference_coefs = np.repeat(1/numref, numref)
    
    if numref == 1:
        ref_part = Model(ref1)
    elif numref == 2:
        ref_part = Model(ref1) + Model(ref2)
    elif numref == 3:
        ref_part = Model(ref1) + Model(ref2) + Model(ref3)
    elif numref > 3:
        ref_part = Model(ref1) + Model(ref2) + Model(ref3)
        print("Warning: Only three reference spectra can currently be fit")
    else:
        ref_part = Model(null_model)


    


    overall_model = ref_part + bg_part + scatter_part + gaussian_part
    params=overall_model.make_params()
    print("Fit parameters")
    print(params)



    if 'bg_val' in params:
        params.add('bg_val', 0.0)
    if 'scatter_int' in params:
        params.add('scatter_int', 1.0, min=0)
        params.add('scatter_delta', 0.5, max=1.0)
        params.add('scatter_power', -4)
        params.add('scatter_a', expr='scatter_delta/(300**scatter_power)') #delta allows inequality that prevents Scatter function returning NaN's
        
        if hold_scatter_power:
          params['scatter_power'].vary = False
        else:
          params['scatter_power'].min = -4
          params['scatter_power'].max = 4

    if 'amp' in params:
        params.add('amp', 0.001)
        params.add('cen', 220.0)
        params.add('wid', 400.0)
    if 'refcoef1' in params:
        params.add('refcoef1', 1.0/(numref+0.001))
        params['refcoef1'].min=0.0
    if 'refcoef2' in params:
        params.add('refcoef2', 1.0/(numref+0.001))
        params['refcoef2'].min=0.0
    if 'refcoef3' in params:
        params.add('refcoef3', 1.0/(numref+0.001))
        params['refcoef3'].min=0.0
   
    print("Parameters set, preparing to fit background.")
    print(mask)
    print(fit_x)
    print(fit_y)
    result = overall_model.fit(fit_y, params, x = fit_x, fit_kws={'maxfev': 200})
    print("Result.report_fit()")
    print(result.fit_report())
    print("report_fit(result)")
    report_fit(result)
    print("report_fit(result.params")
    report_fit(result.params)
    print("report_error(result)")
    report_errors(result)
    #result.plot_fit()
    components = result.eval_components(fit_x=fit_x)
    self.last_fit_result = list([result, components, fit_x, fit_y])
    return(result, components, fit_x, fit_y)

      

        
    


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
class spectral_collection:
  def __init__(self, filename="", load=False):
    self.number_spectra=0
    self.number_files_loaded=0
    self.spectra=list()
    self.filenames=list()
    self.sample_names=list()
    self.rawdata=list()
    self.df = pd.DataFrame()
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
    self.rawdata.append(pd.read_csv(filename, header=None))
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
      xvalues=pd.to_numeric(self.rawdata[self.number_files_loaded-1].loc[2:num_wavelength+1, count*2])
      yvalues=pd.to_numeric(self.rawdata[self.number_files_loaded-1].loc[2:num_wavelength+1, (count*2)+1])
      xvalues.index = np.arange(0, len(xvalues))
      yvalues.index = np.arange(0, len(yvalues))
      
      this_spectrum.wavelength = xvalues
      this_spectrum.absorbance = yvalues
      this_spectrum.sample_name=(self.rawdata[self.number_files_loaded-1].loc[0, (count*2)])
      self.sample_names.append(this_spectrum.sample_name)
      self.spectra.append(this_spectrum)
      self.number_spectra = self.number_spectra+1

      if generatedf:
        xvalues.name="WL:"+names[count*2]
        yvalues.name="AU:"+names[count*2]
        #print(type(xvalues))
        
        tempdf = pd.concat([xvalues,yvalues], axis=1, ignore_index=True)
        tempdf.columns=[xvalues.name, yvalues.name] # not sure why this is necessary, but seems to be
        #print(tempdf)
        if count==0:
          self.df=tempdf
        else:
          self.df = pd.concat([self.df,tempdf], axis=1)
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
