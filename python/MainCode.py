#%%
import cilantro_ancho
from datetime import datetime
from glob import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#input_log = '../data/5Ghz/cilantro_ancho_log5_230519_22_00_07.log'
input_logs = sorted(glob('/Users/vonw/data/icecaps/sleigh/gpr/5Ghz/*.log'))
input_logs = input_logs[105:]

# %% read raw data files
time = []
DM   = []
f    = []
for input_log in input_logs:
    rawdata = cilantro_ancho.logger_parse_RLH(input_log)
    # ....Times in the input_log files are different from the filenames...
    #time.append(gpr.time_DM.values)
    # So, decode time using filename...
    time.append(datetime(2000+int(input_log.split('/')[-1][20:22]), int(input_log.split('/')[-1][22:24]), int(input_log.split('/')[-1][24:26]), int(input_log.split('/')[-1][27:29]), int(input_log.split('/')[-1][30:32]), 0))
    DM.append(rawdata.DM.values)
    # ....Sample delay = 3.867e-09 !! offset of 3.8 ns
    f.append(rawdata.SamplesPerSecond)

gpr = xr.Dataset(data_vars=dict(
                    DM=(['step', 'time'], np.array(DM).T), 
                    f=(['time'], np.array(f))
                ), 
                coords=dict(
                    step=rawdata.step, 
                    time=np.array(time)
                ),
                attrs=rawdata.attrs)

# %% lets plot the data at this stage
gpr.DM.plot()
plt.gca().invert_yaxis()

# %% sampling frequency correction
NumberOfSteps,NumberOfTimes = gpr.DM.shape
tot_time=NumberOfSteps/gpr.f
e1=cilantro_ancho.e_snowdry(350,5e9,-10)    # estimate of dielectric constant
v=3e8/np.sqrt(np.real(e1))                  # velocity in snow
dsnow=v*tot_time/2                          # depth in snow
d0=gpr.SampleDelay*v/2                      # sample delay offset

dscale = np.linspace(d0,dsnow,NumberOfSteps)
max_d  = dscale[-1,:].max()

# ....interpolate to correct for different sampling frequencies
dnew=np.linspace(d0,max_d,NumberOfSteps)

gpr.coords['depth'] = (('depth'), dnew)
gpr['DM2'] = (('depth', 'time'), np.array([np.interp(dnew, dscale[:,time], gpr.DM[:,time]) for time in range(NumberOfTimes)]).T)

# %%
gpr.DM2.plot()
plt.gca().invert_yaxis()

# %%
Ix = np.where(dnew<=3.74)[0]
gpr['DM2'] = gpr.DM2[Ix, :]
dnew = dnew[Ix]

# %% bandpass filter
from statistics import median
from scipy.signal import butter, lfilter, freqz

# ....Copied from https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter
#             and https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html)
###
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
###

dt         = median(np.diff(dnew)*2./v)     # get sampling interval in s
fs         = 1/dt                           # sampling rate in Hz
wavelength = 5.4e9
fmin       = wavelength-3e9
fmax       = wavelength+3e9
order      = 3

gpr['DM3'] = (('depth', 'time'), np.array([butter_bandpass_filter(gpr.DM2[:, time], fmin, fmax, fs, order) for time in range(NumberOfTimes)]).T)

# %%
gpr.DM3.plot()
plt.gca().invert_yaxis()

# %%
