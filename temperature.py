# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

path = r'C:\Users\PSClo\OneDrive\Miscallaneous\Documents\Imperial_physics\Mastersproj\thermo_data'

Tsc_f = np.loadtxt(path + r'\tempscan_f.txt')

Tsc_T = np.loadtxt(path + r'\tempscan_temps.txt')
Tsc_v = np.loadtxt(path + r'\tempscan_v.txt')
Tsc_time = np.loadtxt(path + r'\tempscan_times.txt')

Tsc_time_list = Tsc_time.tolist()
Tsc_v_list = Tsc_v.tolist()
   

for i in range(len(Tsc_time_list)):
    for j in range(len(Tsc_time_list[0])):
        moveindex = len(Tsc_time_list[0])
        movevalue = Tsc_time_list[0][moveindex - 1]
        Tsc_time_list[i][j] = Tsc_time_list[i][j] + i*movevalue
#%%

def flatten(List):
    flatlist = []
    for sublist in List:
        for num in sublist:
            flatlist.append(num)
    return flatlist

flat_timelist = flatten(Tsc_time_list)
flat_vlist = flatten(Tsc_v_list)
flat_freqlist = flatten(Tsc_f)

       
df_f = pd.DataFrame(Tsc_f)
df_v = pd.DataFrame(Tsc_v)
df_times = pd.DataFrame(Tsc_time_list)
#df_f.insert(loc=0, column='temps', value=Tsc_T)


df_list = [df_f, df_v, df_times]
for i, df in enumerate(df_list):
    df_list[i] = df.T
    
df_f = df_list[0]
df_v = df_list[1]
df_times = df_list[2]
N_scans = len(df_list[0].columns)



#%% GLP fitting function for 3 peaks

p0_GLP3 = [192288, 192310, 192330, 
           6,6,6,
           11,11,11,
           0.004, 0.004, 0.004
           ]


def GLP3(w, 
        c1, c2, c3, 
        deltaw1, deltaw2, deltaw3, 
        A1, A2, A3,
        m1, m2, m3,
        D=0):
    """
    w: (array) frequencies
    c: (float) centre frequency
    deltaw: (float) full width at half maximum
    A: (float) amplitude/height of peak
    m: (float) mixing ratio; m=0 --> pure Gaussian, m=1 --> pure Lorentzian
    D: (float, optional) vertical displacement
    """
    
    GLP1 = (A1*np.exp(-4*np.log(2)*(1-m1)*(w-c1)**2/deltaw1**2) * 
            1/(1+4*m1*(w-c1)**2/deltaw1**2))
    GLP2 = (A2*np.exp(-4*np.log(2)*(1-m2)*(w-c2)**2/deltaw2**2) * 
            1/(1+4*m2*(w-c2)**2/deltaw2**2))
    GLP3 = (A3*np.exp(-4*np.log(2)*(1-m3)*(w-c3)**2/deltaw3**2) * 
            1/(1+4*m3*(w-c3)**2/deltaw3**2))

    return GLP1 + GLP2 + GLP3


def fitting_GLP3(freqs, voltages):
    
    GLP_fit_3 = curve_fit(GLP3, freqs, voltages, p0 = p0_GLP3, maxfev = 50000)
    p_GLP3_err = np.sqrt(np.diag(GLP_fit_3[1]))
    Varray_fitGLP3 = GLP3(freqs, *GLP_fit_3[0])
    peak_Varray = [GLP_fit_3[0][6], GLP_fit_3[0][7], GLP_fit_3[0][8]]
    central_Freqs = [GLP_fit_3[0][0], GLP_fit_3[0][1], GLP_fit_3[0][2]]
    fwhm_Vals = [GLP_fit_3[0][3], GLP_fit_3[0][4], GLP_fit_3[0][5]]
    m_vals = [GLP_fit_3[0][9], GLP_fit_3[0][10], GLP_fit_3[0][11]]
    
    
    return(GLP_fit_3, peak_Varray, central_Freqs, fwhm_Vals, m_vals)


def GLP4(w, 
        c1, c2, c3, c4,
        deltaw1, deltaw2, deltaw3, deltaw4 ,
        A1, A2, A3,A4,
        m1, m2, m3, m4,
        D=0):
    """
    w: (array) frequencies
    c: (float) centre frequency
    deltaw: (float) full width at half maximum
    A: (float) amplitude/height of peak
    m: (float) mixing ratio; m=0 --> pure Gaussian, m=1 --> pure Lorentzian
    D: (float, optional) vertical displacement
    """
    
    GLP1 = (A1*np.exp(-4*np.log(2)*(1-m1)*(w-c1)**2/deltaw1**2) * 
            1/(1+4*m1*(w-c1)**2/deltaw1**2))
    GLP2 = (A2*np.exp(-4*np.log(2)*(1-m2)*(w-c2)**2/deltaw2**2) * 
            1/(1+4*m2*(w-c2)**2/deltaw2**2))
    GLP3 = (A3*np.exp(-4*np.log(2)*(1-m3)*(w-c3)**2/deltaw3**2) * 
            1/(1+4*m3*(w-c3)**2/deltaw3**2))
    GLP4 = (A4*np.exp(-4*np.log(2)*(1-m4)*(w-c4)**2/deltaw4**2) * 
            1/(1+4*m4*(w-c4)**2/deltaw4**2))

    return GLP1 + GLP2 + GLP3 + GLP4


def fitting_GLP4(freqs, voltages):
    
    GLP_fit_4 = curve_fit(GLP4, freqs, voltages, p0 = p0_GLP4, maxfev = 500000)
    p_GLP4_err = np.sqrt(np.diag(GLP_fit_4[1]))
    Varray_fitGLP3 = GLP4(freqs, *GLP_fit_4[0])
    peak_Varray = [GLP_fit_4[0][6], GLP_fit_4[0][7], GLP_fit_4[0][8]]
    central_Freqs = [GLP_fit_4[0][0], GLP_fit_4[0][1], GLP_fit_4[0][2]]
    fwhm_Vals = [GLP_fit_4[0][3], GLP_fit_4[0][4], GLP_fit_4[0][5]]
    m_vals = [GLP_fit_4[0][9], GLP_fit_4[0][10], GLP_fit_4[0][11]]
    
    
    return(GLP_fit_4, peak_Varray, central_Freqs, fwhm_Vals, m_vals)



def GLP2(w, 
        c1, c2, 
        deltaw1, deltaw2,
        A1, A2, 
        m1, m2, 
        D=0):
    """
    w: (array) frequencies
    c: (float) centre frequency
    deltaw: (float) full width at half maximum
    A: (float) amplitude/height of peak
    m: (float) mixing ratio; m=0 --> pure Gaussian, m=1 --> pure Lorentzian
    D: (float, optional) vertical displacement
    """
    
    GLP1 = (A1*np.exp(-4*np.log(2)*(1-m1)*(w-c1)**2/deltaw1**2) * 
            1/(1+4*m1*(w-c1)**2/deltaw1**2))
    GLP2 = (A2*np.exp(-4*np.log(2)*(1-m2)*(w-c2)**2/deltaw2**2) * 
            1/(1+4*m2*(w-c2)**2/deltaw2**2))
    

    return GLP1 + GLP2 


def fitting_GLP2(freqs, voltages):
    
    GLP_fit_2 = curve_fit(GLP2, freqs, voltages, p0 = p0_GLP2)
    p_GLP2_err = np.sqrt(np.diag(GLP_fit_2[1]))
    Varray_fitGLP2 = GLP2(freqs, *GLP_fit_2[0])
    peak_Varray = [GLP_fit_2[0][6], GLP_fit_2[0][7]]
    central_Freqs = [GLP_fit_2[0][0], GLP_fit_2[0][1]]
    fwhm_Vals = [GLP_fit_2[0][3], GLP_fit_2[0][4]]
    
    
    return(GLP_fit_2, peak_Varray, central_Freqs, fwhm_Vals)




 #%% Fitting Data 

voltage_peaks = []
central_freqs = []
fwhm_vals = []

for i in range(N_scans):
    x = fitting_GLP3(df_f[i], df_v[i])
    voltage_peaks.append(x[1])
    central_freqs.append(x[2])
    fwhm_vals.append(x[3])
   

#%% 


mean_voltage_peaks = []
FSRs = []
mean_FSRs = []
mean_fwhms = []
for i in range(N_scans):
    FSR = []
    for j in range(len(central_freqs[0]) -1):
        fsr = central_freqs[i][j+1] - central_freqs[i][j]
        FSR.append(fsr)
    FSRs.append(FSR)
    
    mean_fsr = np.mean(FSRs[i])
    mean_FSRs.append(mean_fsr)
    
    mx = np.mean(voltage_peaks[i])
    mean_voltage_peaks.append(mx)
    
    mean_fwhm = np.mean(fwhm_vals[i])
    mean_fwhms.append(mean_fwhm)
    
flat_FSR = flatten(FSRs)
mean_FSRs = []
for i in range(len(FSRs)):
    mean_fsr = np.mean(FSRs[i])
    mean_FSRs.append(mean_fsr)


my_data_path = r'C:\Users\PSClo\OneDrive\Miscallaneous\Documents\Imperial_physics\Mastersproj\thermo_data' 

#Save both CSV and txt file
df_labtemp_data = pd.DataFrame(list(zip(mean_voltage_peaks, mean_FSRs, mean_fwhms)))
df_labtemp_data.columns = ['mean_voltage_peaks', 'mean_FSRs', 'mean_FWHM']
df_labtemp_data.to_csv(path + '\df_lab_temp_data.csv', index = False)
np.savetxt(path + '\df_lab_temp_data.txt', df_labtemp_data, delimiter = ',')

#%% Plotting Data
x_A = np.loadtxt(path + '\df_lab_temp_data.csv', delimiter = ',', skiprows =1)
temp_data = pd.DataFrame(x_A, columns = ['mean_voltage_peaks', 'mean_FSRs', 'mean_FWHM'])



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

plt.plot(Tsc_T, temp_data['mean_FSRs'], '-o', c='blue', mfc = 'yellow', mec = 'k')
plt.xlabel("Temperature (K)")
plt.ylabel("FSR")
plt.title('Temperature measurement - Temperature v.s FSR')
plt.show()

plt.plot(Tsc_T, temp_data['mean_voltage_peaks'],'-o', c='blue', mfc = 'red', mec = 'k')
plt.xlabel("Temperature (K)")
plt.ylabel("Voltage (V)")
plt.title('Temperature measurement - Temperature v.s Peak Voltage')
#plt.text(0.5,0.5, "Variance = {}".format(np.var(mean_voltage_peaks), "%.3f"))
plt.show()

#%%
first_cfreq = 192284
freq_shift = -0.525
p = []
voltage_peaks = []
central_freqs = []
fwhm_vals = []
m_vals = []
for i in range(len([5,5])):
    
    p0_GLP4 = [first_cfreq + i*freq_shift -1, first_cfreq + 20 + i*freq_shift - 1, first_cfreq + 40 + i*freq_shift -1 , first_cfreq + 60 + i*freq_shift -1,
               6, 6, 6, 6, 
               12, 12, 12, 12,
               0.004, 0.004, 0.004, 0.004
               ]
    p.append(p0_GLP4)
    
    x = fitting_GLP3(freq_data[i], volt_data[i])
    voltage_peaks.append(x[1])
    central_freqs.append(x[2])
    fwhm_vals.append(x[3])
    m_vals.append(x[4])
    
    
    
    
    













