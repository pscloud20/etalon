# -*- coding: utf-8 -*-
"""
18 Oct 2022
Phil, Pablo, and Steve

Here is some initial code to try a computer-controlled transmission experiment.
"""

import sys # Used for specifying system paths to other modules
newPaths = ('G:\\Code', 'G:\\Code\\Instruments','G:\\Code\\Instruments\\Wavemeter','G:\\Code\\Experiments\\Experiment_Tools',
            'G:\\Code\\Instruments\\PicoScope', 'G:\\Code\\Instruments\\Swabian', 'G:\\Code\\Development\\20220823 SL TR auto',
            'G:\\Code\\Development\\20220825 SL TR multi-peak')
for i in newPaths:
   if i not in sys.path: sys.path.append(i) 
   

from sklab_instruments import wavemeter, EXFO, nuPhotonEDFA, labJackT4

from time import sleep

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

########

"""
wm = wavemeter()
answer = wm.readF(1)
print(answer)
wm.close()
"""

"""
exfo = EXFO(p=32)
answer = exfo.readI()
print(answer)
exfo.close()
"""

"""
lj = labJackT4('LabJack1')
answer = lj.readV(3)
print(answer)

lj.close()
"""

def resonanceScan(freqs, deltat, #single scan params.
                  n_scans=1, deltat_scans=0, #between sans params,
                  savedata=False):
    """
    1.) Creates instruments and sets laser frequency for each frequency value in freqs list 
    2.) Creates empty lists and loops through frequencies and reads the voltage at that frequency
        - like the pause between scans we need a pause between measurements for the frequency to be set (lower limit ~ 0.8)
    3.) Appends values to empty lists and then closes all instruments
    
    
    freqs : list (this allows us to exactly specify what frequencies we would like to input)
    deltat: array((float)) time step between data points. Here sleeps laser so that
            deltaf is made and frequency stabalizes.
    deltat_scans : length 
    n_scans: (int): number of (identical) scans to take.
    deltat_scans: (float) time step between scans minus min time between scans
                    which is (n+1)*sleep_t by construcion.
    savedata: (True/False) whether to savedata as .txt or not.
    
    """

        
    voltages2 = []
    frequencies2 = []
    
    for i in range(len(freqs)):
        #create instruments - wavemeter, laser, labjack
        wm = wavemeter()
        exfo = EXFO(p=32)
        lj = labJackT4('LabJack1')
              
        exfo.setF(freqs[i][0])
        sleep(deltat_scans) #pause required such that there is sufficient time for the laser to change frequency output

        voltages1 = []
        frequencies1 = []
        for j in freqs[i]:
            exfo.setF(j)
            sleep(deltat[i])
            f = wm.readF(1)
            v = lj.readV(3)
            voltages1.extend([v])
            frequencies1.extend([f])
            
        #close instruments
        wm.close()
        exfo.close()
        lj.close()
        
        voltages2.extend([voltages1])
        frequencies2.extend([frequencies1])
        
            
        # if savedata == True:
            
        #     np.savetxt(r'G:\Code\Development\20221018 intro device use\Data\frequencies_data_p{}.txt'.format(laser_current), frequencies, header = 'n = {}, $\Delta$f = {}, freq_init = {}, laser_current = {}'.format(n, deltaf, freq_i, laser_current))
        #     np.savetxt(r'G:\Code\Development\20221018 intro device use\Data\voltages_data_p{}.txt'.format(laser_current), voltages)
        # else:
        #     pass
    
    return(frequencies2, voltages2) 


# =============================================================================

def stabilityScan(deltat, n, freq_fix, #between points params.
                  deltat_scans=1 #between sans params.
                  ): 
    """
    1.) Gives the laser a set frequency 
    2.) Creates a time axis list 
    3.) Measures voltage at each time point in the time axis list
    
    
    deltaf: (array(float)) frequency step size.
    n: (array(int)) number of (time) steps to take (so intervals, not points)
    freq_fix : (array(float)) fixed frequency which we are measuring at 
    deltat: array((float)) time step between data points. Here sleeps laser so that
            deltaf is made and frequency stabalizes.
    n_scans: (int): number of (identical) scans to take.
    deltat_scans: (float) time step between scans minus min time between scans
                    which is (n+1)*sleep_t by construcion.
    """
    
    voltages2 = []
    time2 = []
    
    for i in range(len(freq_fix)):
        
        wm = wavemeter()
        exfo = EXFO(p=32)
        lj = labJackT4('LabJack1')
          
        time1 = np.arange(0, n[i]*deltat[i] + deltat[i], deltat[i])
        time2.extend([time1]) # might be better to append arrays not lists?
                
        voltages1 = []
        exfo.setF(freq_fix[i])
        sleep(deltat_scans)
        for j in time1:
            sleep(deltat[i])
            v = lj.readV(3)
            voltages1.extend([v])
        voltages2.extend([voltages1])
        
        wm.close()
        exfo.close()
        lj.close()
        

    return(time2, voltages2)
# =============================================================================
        

def plot_resonanceScan(printdatainfo):
    
    #fdata = np.loadtxt(r'G:\Code\Development\20221018 intro device use\Data\frequencies_data.txt')
    #vdata = np.loadtxt(r'G:\Code\Development\20221018 intro device use\Data\voltages_data.txt')
    #fdata1 = pd.DataFrame(fdata)
    fdata = np.loadtxt(r'G:\Code\Development\20221018 intro device use\Data\frequencies_data_p{}.txt'.format(laser_current))
    vdata = np.loadtxt(r'G:\Code\Development\20221018 intro device use\Data\voltages_data_p{}.txt'.format(laser_current))
    
    
    plt.rcParams.update(params)
    
    plt.plot(np.array(fdata), vdata, marker = 'x')
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Voltage (mV)")
    plt.title('Plot of Resonance Peaks for a laser in F-P cavity')
    plt.grid()
    #plt.xscale('linear')
    
    
    if printdatainfo == True:
        
        plt.figtext(x= 0, y =0 , s = 'n = {}, $\Delta$f = {}, freq_init = {}, laser_current = {}'.format(n, deltaf, freq_i, laser_current), fontsize = 10)
    else:
        pass
    return


#%%
# plot_resonanceScan(printdatainfo=True)

# deltaf=0.5
# n=150
# freq_i = 192280
# laser_current  = 40
# sleep_t = 1

# datagen = resonanceScan(deltaf, n, freq_i, savedata = True, sleep_t = sleep_t)
# #%%
# plt.scatter(datagen[1], datagen[0])
# #%%
# stabdatagen = stabilityScan(delta_t = 10, N = 5, freq_fix = 192295.724 , sleep_t = 1)



# =============================================================================
# #%%
# 
# filef = open(r'G:\Code\Development\20221018 intro device use\Data\frequencies_data.txt', 'r')
# fdata1 = filef.readline()
# 
# #%%
#     
# 
# resonanceScan(deltaf = 0.5, n = 155, freq_i = 192280, savedata = False)
# 
# fdata = np.loadtxt(os.getcwd() + r'\Data\frequencies_data.txt')
# vdata = np.loadtxt(os.getcwd() + r'\Data\voltages_data.txt')
# 
# plot_resonanceScan(fdata, vdata)
# #%%
# 
# wm = wavemeter()
# exfo = EXFO(p=32)
# lj = labJackT4('LabJack1')
# 
# deltaf=0.5
# n=155
# freq_i = 192280
# freqs = np.arange(freq_i, freq_i+n*deltaf, deltaf)
# 
# voltages = []
# frequencies = []
# for i in freqs:
#     exfo.setF(i)
#     sleep(1)
#     f = wm.readF(1)
#     # print(f)
#     v = lj.readV(3)
#     voltages.append(v)
#     frequencies.append(f)
#     # print(v)
# 
# wm.close()
# exfo.close()
# lj.close()
# 
# #%%
# plt.plot(np.array(frequencies), voltages, marker = "x")
# plt.xlabel("Frquency (GHz)")
# plt.ylabel("Voltage (mV)")
# plt.title('Plot of Resonance Peaks for a laser in F-P cavity')
# plt.grid()
# 
# np.savetxt('frequencies_data.txt', frequencies)
# np.savetxt('voltages_data.txt', voltages)
# 
# plt.xticks(np.arange(freq_i,  freq_i+n*deltaf, 10))
# plt.tight_layout()
# 
# #%% Analysis
# 
# frequencies_peak1 = [i for i in frequencies if i < 192290]
# voltages_peak1 = voltages[:len(frequencies_peak1)]
# max_peak1 = max(voltages_peak1)
# freq_peak1 = frequencies_peak1[voltages_peak1.index(max_peak1)]
# print(max_peak1, freq_peak1)
# plt.plot(frequencies_peak1, voltages_peak1, 'x')
# plt.axvline(freq_peak1, color = 'red')
# 
# 
# 
# #%%
# 
# wm = wavemeter()
# exfo = EXFO(p=32)
# lj = labJackT4('LabJack1')
# 
# exfo.setF(192303)
# f = wm.readF(1)
# print(f)
# 
# wm.close()
# exfo.close()
# lj.close()
# #%%
# wm = wavemeter()
# exfo = EXFO(p=32)
# 
# frequencies = []
# for i in freqs:
#     exfo.setF(i)
#     sleep(3)
#     f = wm.readF(1)
#     print(f)
#     frequencies.append(f)
#     #print(answer)
# 
# wm.close()
# exfo.close()
# #%%
# plt.plot(np.arange(0, n, 1), frequencies)
# 
# =============================================================================
