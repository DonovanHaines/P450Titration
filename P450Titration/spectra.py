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
    

  def baseline_correct(self, reference_spectra=None, control=None):
      default_control = { #default values
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
                "wavelength_low": self.wl_low,
                "wavelength_high": self.wl_high,
                "run_all": False}
      if control==None:
          control = default_control
      else:
          #TODO check for missing vital paramters here
          if "wavelength_low" not in control:
              control["wavelength_low"] = default_control["wavelength_low"]
          if "wavelength_high" not in control:
              control["wavelength_high"] = default_control["wavelength_high"]
      
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
      wave_low = control.get("wavelength_low", 0.0)
      wave_high = control.get("wavelength_high", 1500.0)
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
  
      def gaussian(x, gaussian_amp, gaussian_cen, gaussian_wid):
        #"1-d gaussian: gaussian(x, amp, cen, wid)"
        return np.asarray([gaussian_amp * np.exp(-(xt-gaussian_cen)**2.0 / gaussian_wid) for xt in x])
     
      def null_model(x):
          return np.zeros(len(x))


      if control["fit_gaussian"]: 
          gaussian_part = Model(gaussian)
      else: 
          gaussian_part = Model(null_model)

      if control["fit_scatter"]:
          scatter_part = Model(scatter)
      else:
          scatter_part = Model(null_model)
    
      if control["fit_bg"]:
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

      # now to initialize all the parameters
      quick_not = lambda item: not item

      if 'bg_val' in params:
          params['bg_val'].set(value = control.get("bg_val"))
          params['bg_val'].set(vary = quick_not(control.get("bg_val_fixed", False)))        
          params['bg_val'].set(min = control.get("bg_val_min", -math.inf))
          params['bg_val'].set(max = control.get("bg_val_max", +math.inf))

      if 'scatter_int' in params:
          params['scatter_int'].set(value=control.get("scatter_int"), vary=(not control.get('sscatter_int', False)))
          params['scatter_int'].set(vary = quick_not(control.get("scatter_int_fixed", False)))  
          params["scatter_int"].set(min=control.get("scatter_int_min", -math.inf))
          params["scatter_int"].set(max=control.get("scatter_int_max", +math.inf))

      if 'scatter_a' in params:
          params['scatter_a'].set(value=control.get("scatter_a"))
          params['scatter_a'].set(vary = quick_not(control.get("scatter_a_fixed", False)))  
          params["scatter_a"].set(min=control.get("scatter_a_min", -math.inf))
          params["scatter_a"].set(max=control.get("scatter_a_max", +math.inf))
          #  and add scatter_delta which enables to use of the expression to limit scatter_a to meaninful values
          params.add("scatter_delta", value=control.get("scatter_delta"))
          params["scatter_delta"].set(min=control.get("scatter_delta_min", -math.inf))
          if "scatter_delta_max" in control:
              params["scatter_delta"].set(max=control.get("scatter_delta_max", +math.inf))
          else:
              params["scatter_delta"].set(max=1.0)
      if 'scatter_power' in params:
          params['scatter_power'].set(value=control.get("scatter_power"))
          params['scatter_power'].set(vary = quick_not(control.get("scatter_power", False)))  
          params["scatter_power"].set(min=control.get("scatter_power_min", -math.inf))
          params["scatter_power"].set(max=control.get("scatter_power_max", +math.inf))
          params.add('scatter_a', expr='scatter_delta/(300**scatter_power)') #delta allows inequality that prevents Scatter function returning NaN's


      if 'gaussian_amp' in params:
          params['gaussian_amp'].set(value=control.get("gaussian_amp"))
          params['gaussian_amp'].set(vary = quick_not(control.get("gaussian_amp_fixed", False)))
          params["gaussian_amp"].set(min=control.get("gaussian_amp_min", -math.inf))
          params["gaussian_amp"].set(max=control.get("gaussian_amp_max", +math.inf))
      if 'gaussian_wid' in params:
          params['gaussian_wid'].set(value=control.get("gaussian_wid"))
          params['gaussian_wid'].set(vary = quick_not(control.get("gaussian_wid_fixed", False)))
          params["gaussian_wid"].set(min=control.get("gaussian_wid_min", -math.inf))
          params["gaussian_wid"].set(max=control.get("gaussian_wid_max", +math.inf))
      if 'gaussian_cen' in params:
          params['gaussian_cen'].set(value=control.get("gaussian_cen"))
          params['gaussian_cen'].set(vary = quick_not(control.get("gaussian_cen_fixed", False)))
          params["gaussian_cen"].set(min=control.get("gaussian_cen_min", -math.inf))
          params["gaussian_cen"].set(max=control.get("gaussian_cen_max", +math.inf))

      if 'refcoef1' in params:
          params['refcoef1'].set(value=control.get("ref_coef1"))
          params["refcoef1"].set(min=control.get("ref_coef1_min", -math.inf))
          params["refcoef1"].set(max=control.get("ref_coef1_max", +math.inf))

      if 'refcoef2' in params:
          params['refcoef2'].set(value=control.get("ref_coef2"))
          params["refcoef2"].set(min=control.get("ref_coef2_min", -math.inf))
          params["refcoef2"].set(max=control.get("ref_coef2_max", +math.inf))

      if 'refcoef3' in params:
          params['refcoef3'].set(value=control.get("ref_coef3"))
          params["refcoef3"].set(min=control.get("ref_coef3_min", -math.inf))
          params["refcoef3"].set(max=control.get("ref_coef3_max", +math.inf))
   
      print("Parameters set, preparing to fit background.")
      print(mask)
      print(fit_x)
      print(fit_y)
      result = overall_model.fit(fit_y, params, x = fit_x, fit_kws={'maxfev': 200}, iter_cb=self.iter_monitor)
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
      self.last_fit_result = list([result, components, fit_x, fit_y, control])
      
      return (result, components, fit_x, fit_y, control)
  
  ######################################################################
  ######################################################################
  def iter_monitor(self, params, iter, resid, *args, **kws):
      print("Iteration no. {0} SSR: {1}".format(iter, sum(resid*resid)))
      

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
      this_spectrum.wl_low = min(xvalues)
      this_spectrum.wl_high = max(xvalues)
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
