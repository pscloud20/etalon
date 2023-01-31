# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:29:46 2022

@author: PSClo
"""
import glob
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

path = r'C:\Users\PSClo\OneDrive\Miscallaneous\Documents\Imperial_physics\Mastersproj\thermo_data'

T = np.loadtxt(path + r'\multi_scans\heatcan_temps299.600.txt')
v = np.loadtxt(path + r'\multi_scans\heatscan_v299.600.txt')
f = np.loadtxt(path + r'\multi_scans\heatscan_f299.600.txt')

files = []
temp_files = glob.glob(path + r'\multi_scans\heatcan_temps*.txt')
freq_files = glob.glob(path + r'\multi_scans\heatscan_f*.txt')
volt_files = glob.glob(path + r'\multi_scans\heatscan_v*.txt') 
time_files = glob.glob(path + r'\multi_scans\heatscan_times*.txt')
    
def dataimport_tolist(T_files, f_files, v_files, t_files):
    dataA = []
    dataB = []
    dataC = []
    dataD = []
    T_data = []
    freq_data = []
    volt_data = []
    time_data = []
    for f in range(len(temp_files)):
        T = np.loadtxt(T_files[f])
        F = np.loadtxt(f_files[f])
        V = np.loadtxt(v_files[f])
        t = np.loadtxt(t_files[f])
        dataA.append(T)
        dataB.append(F)
        dataC.append(V)
        dataD.append(t)
        T = dataA[f].tolist()
        F = dataB[f].tolist() 
        V = dataC[f].tolist()
        T_data.append(T)
        freq_data.append(F)
        volt_data.append(V)
        time_data.append(t)
    
    if (np.shape(volt_data) == np.shape(freq_data) == np.shape(time_data)) == True:
        print('import successful')
    
    return(T_data, freq_data, volt_data, time_data)

T_data, freq_data, volt_data, time_data = dataimport_tolist(temp_files, freq_files, volt_files, time_files)


#%%
first_cfreq = 192284
freq_shift = -0.525
p = []
voltage_peaks = []
central_freqs = []
fwhm_vals = []
m_vals = []
for i in range(len(T_data)):
    
    p0_GLP3 = [first_cfreq + i*freq_shift -1, first_cfreq + 20 + i*freq_shift - 1, first_cfreq + 40 + i*freq_shift -1,
               6, 6, 6,
               12, 12, 12, 
               0.004, 0.004, 0.004
               ]
    p.append(p0_GLP3)
    
    x = fitting_GLP3(freq_data[i], volt_data[i])
    voltage_peaks.append(x[1])
    central_freqs.append(x[2])
    fwhm_vals.append(x[3])
    m_vals.append(x[4])
    
plt.scatter(freq_data[1], volt_data[1])
plt.scatter(central_freqs[1], voltage_peaks[1])
    


#%%
p0_GLP3 = [first_cfreq, first_cfreq + 20, first_cfreq,
               6, 6, 6,
               12, 12, 12, 
               0.004, 0.004, 0.004
               ]
p.append(p0_GLP3)
    



x2 = fitting_GLP3(freq_data[2], volt_data[2])
voltage_peaks.append(x[1])
central_freqs.append(x[2])
fwhm_vals.append(x[3])
m_vals.append(x[4])

plt.scatter(freq_data[2], volt_data[2])
plt.scatter(x2[2], x2[1])


    

#%%
params = {
    'font.family' : 'serif',
    'font.size': 13.5,
    'figure.figsize':  [8.8, 8.8/1.618],
    'axes.grid': True,
    'legend.loc' :'best',
    'legend.fontsize': 10,
    'legend.handlelength':2,

}

plt.rcParams.update(params)

plt.plot(time_data[1], volt_data[1], '-o', c='blue', mfc = 'yellow', mec = 'k')
plt.show()

#%%
n = 12


   
p0_GLP3 = [192280, 192300, 192320, 
           6,6,6,
           12,12,12,
           0.004, 0.004, 0.004
           ]

def fit_lists_temp_analysis(freq_data, volt_data):
    voltage_peaks = []
    central_freqs = []
    fwhm_vals = []
    m_vals = []
    
    for i in range(n):
        x = fitting_GLP3(freq_data[i], volt_data[i])
        voltage_peaks.append(x[1])
        central_freqs.append(x[2])
        fwhm_vals.append(x[3])
        m_vals.append(x[4])
        
    fit_vals = []
    fit_freqs = np.arange(192279, 192340, 0.0001)
    
    
        
        
    return(voltage_peaks, central_freqs, fwhm_vals, m_vals)
        
voltage_peaks, central_freqs, fwhm_vals, m_vals = fit_lists_temp_analysis(freq_data, volt_data)



def plot_temp_analysis_fit(plot, subplot, fitplot):
    
    

    if plot == True:
        colors = plt.cm.cividis(np.linspace(0,1, n))
    
        for i in range(n):
                #plt.subplot(3,4, i+1)
                plt.plot(freq_data[i], volt_data[i], '-', color = colors[i])
                plt.title('Resonance peak shifting due to temp change')
                plt.xlabel('Frequency (GHz)')
                plt.ylabel('Voltage (V)')
                plt.plot(central_freqs[i], voltage_peaks[i], 'x', color = 'red')
        plt.show()
                
    if subplot == True:
        for i in range(n):
            colors = plt.cm.cividis(np.linspace(0,1, n))
            plt.subplot(3,4, i+1)
            plt.plot(freq_data[i], volt_data[i], '-', color = colors[i])
            plt.title('Resonance peak shifting due to temp change')
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Voltage (V)')
            #plt.plot(central_freqs[i], voltage_peaks[i], 'x', color = 'red')
        plt.show()
    
    if fitplot == True:
        fit_vals = []
        fit_freqs = np.arange(192279, 192340, 0.0001)
        
        for i in range(n):
        
            fitting = GLP3(freq_data[i], 
                           central_freqs[i][0], central_freqs[i][1], central_freqs[i][2],
                           fwhm_vals[i][0], fwhm_vals[i][1], fwhm_vals[i][2],
                           voltage_peaks[i][0], voltage_peaks[i][1], voltage_peaks[i][2],
                           m_vals[i][0], m_vals[i][1], m_vals[i][2])
            fit_vals.append(fitting)
            
        for i in range(n):
            #plt.subplot(3,4, i+1)
            colors = plt.cm.cividis(np.linspace(0,1, n))
            plt.plot(freq_data[i], fit_vals[i], color = colors[i])
            plt.plot(freq_data[i], volt_data[i], 'x', markersize = 3)
        plt.show()
                
    return
    
plot_temp_analysis_fit(plot = True, subplot = False, fitplot = True)
#%%
voltage_peaks = []
central_freqs = []
fwhm_vals = []
m_vals = []

for i in range(n):
    x = fitting_GLP3(freq_data[i], volt_data[i])
    voltage_peaks.append(x[1])
    central_freqs.append(x[2])
    fwhm_vals.append(x[3])
    m_vals.append(x[4])

flat_vpeaks = flatten(voltage_peaks)    
flat_cfreqs = flatten(central_freqs)
    

colors = plt.cm.cividis(np.linspace(0,1, n))

for i in range(n):
    #plt.subplot(3,4, i+1)
    plt.plot(freq_data[i], volt_data[i], '-', color = colors[i])
    plt.title('Resonance peak shifting due to temp change')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Voltage (V)')
    plt.plot(central_freqs[i], voltage_peaks[i], 'x', color = 'red')
    

#plt.plot(central_freqs, voltage_peaks, 'x', color = 'red')
#plt.show()
#%%
fit_vals = []
fit_freqs = np.arange(192279, 192340, 0.0001)

for i in range(n):

    fitting = GLP3(freq_data[i], 
                   central_freqs[i][0], central_freqs[i][1], central_freqs[i][2],
                   fwhm_vals[i][0], fwhm_vals[i][1], fwhm_vals[i][2],
                   voltage_peaks[i][0], voltage_peaks[i][1], voltage_peaks[i][2],
                   m_vals[i][0], m_vals[i][1], m_vals[i][2])
    fit_vals.append(fitting)

for i in range(n):
    #plt.subplot(3,4, i+1)
    plt.plot(freq_data[i], fit_vals[i])
    plt.plot(freq_data[i], volt_data[i], 'x', markersize = 3)
    #plt.show()
 
#%% Finer fit of frequencies
 
fine_fit_vals = []
for i, j in zip(range(n), range(len(fit_freqs))):

    fitting = GLP3(fit_freqs, 
                   central_freqs[i][0], central_freqs[i][1], central_freqs[i][2],
                   fwhm_vals[i][0], fwhm_vals[i][1], fwhm_vals[i][2],
                   voltage_peaks[i][0], voltage_peaks[i][1], voltage_peaks[i][2],
                   m_vals[i][0], m_vals[i][1], m_vals[i][2])
    fine_fit_vals.append(fitting)
    
for i, j in zip(range(len(fit_freqs)), range(n)):
    #plt.subplot(3,4, i+1)
    plt.plot(fit_freqs, fine_fit_vals[i])
    plt.plot(freq_data[j], volt_data[j], 'x', markersize = 3)
    plt.show()

#%%

central = []   
for i in(range(len(central_freqs))):
    central.append(central_freqs[i])
plt.plot(T_data, central, 'x')
               
               
               
    
    
    
    
    
    
    
    
    
    
    
    
    
    
















