import pandas
from tkinter import *

import numpy as np
import math
from scipy.optimize import optimize
from lmfit import Model

import bisect
import random


def gaussian(x, amp, cen, wid):
  #"1-d gaussian: gaussian(x, amp, cen, wid)"
  return amp * np.exp(-(x-cen)**2 / wid)
gaussian_part = Model(gaussian)    
overall_model = gaussian_part
params=overall_model.make_params()
params['amp'].value = 1
params['cen'].value = 418
params['wid'].value = 100

x = np.linspace (300,800,500)
y=overall_model.eval(params, x=x)
result = overall_model.fit(y, params, x = x)
