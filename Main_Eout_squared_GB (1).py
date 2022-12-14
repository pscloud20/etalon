# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 21:50:35 2022

@author: pgara
"""
import numpy as np
import matplotlib.pyplot as plt
from Eout_squared import EGout2, EGout2_2, line, finesse
# from Integration_methods import RK4, Trap, ExtTrap, Euler, AB2
from scipy.signal import find_peaks
from scipy.optimize import curve_fit



import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
#%%
lamda = 1560e-9
k = 2*np.pi/lamda
L = 5e-3
r = np.sqrt(0.75)
d = 2e-3
n_g = 1.5

freqi = 192280
freqf = freqi+20*6 + 3
deltaf =  1e-3# GHz
freq = np.arange(freqi, freqf+deltaf, deltaf) 
rows=4
cols=4

Larray = np.linspace(L, L+1e-6, rows*cols, endpoint = True)
karray = 2*np.pi/3e8 * freq*1e9



#%%


# =============================================================================
# GAUSSIAN BEAM
# =============================================================================



#%%
# =============================================================================
# VARY L, FIX RHO = 0, FIX THETA = 0, WAIST = 1E-3
# =============================================================================

m = np.arange(0, int(69), 1)
Larray3 = np.linspace(L, L+2.5*lamda, int(1e4))
EGout2_L=[]
for i in range(len(Larray3)):
    print("##########", i)
    EGout2i = EGout2(0,0,0,
                L = Larray3[i], 
                d = d,
                n_g = n_g,
                k = k, 
                thetai = 0*np.pi/180, 
                r = r, 
                m = m, 
                w0 = 1e-3, 
                )
    EGout2_L.append(EGout2i)
#%%
plt.figure()
plt.plot(Larray3*1e3, EGout2_L, color = "RED", zorder = len(Larray3))

plt.xlabel("Cavity length (mm)")
plt.ylabel("Amplitude")
# plt.xlim(((L+1.6e-7-6e-8)*1e3, (L+2.2e-7)*1e3))
plt.grid(zorder = 0)
plt.tight_layout()
# plt.savefig("Single_peak_Gaussian_m70_L7.3e-8_L2.2e-7.jpg", dpi = 400)
#%%
# =============================================================================
# VARY L, FIX RHO = 0, FIX THETA = 0, SCAN FREQ, WAIST = 1E-3
# =============================================================================

m = np.arange(0, int(69), 1)
EGout2_L_k=[]
for i in range(len(Larray)):
    print("##########", i)
    EGout2i = EGout2_2(0,0,0, 
                L = Larray[i], 
                d = d,
                n_g = n_g,
                k = karray, 
                thetai = 0, 
                r = r, 
                m = m, 
                w0 = 1e-3, 
                )
    EGout2_L_k.append(EGout2i)
#%%
fig, axs = plt.subplots(rows, cols, sharey=False, sharex = False, figsize=(20, 16))

counter = 0
for j in range(cols):
    
    for i in range(rows):
        
        axs[i, j].plot(freq, EGout2_L_k[counter], color = "r", zorder = len(freq))
        
        axs[rows-1, j].set_xlabel("Frequency (GHz)")
        axs[i, 0].set_ylabel("Amplitude")
        axs[i,j].grid(zorder = 0, alpha = 0.3)
        axs[i,j].set_title(r"$\theta_i$=0°, L={} mm".format(round(Larray[counter]*1e3, 5)))
        
        fig.tight_layout(h_pad = 2)
        
        counter += 1
        
# fig.savefig("Gaussian_varyL_fixRho0_fixTheta0,fixWaist1e-3.jpg", dpi = 400)

#%%
# =============================================================================
# FIX L, FIX RHO = 0, VARY THETA, SCAN FREQ, WAIST = 1E-3
# =============================================================================

from Eout_squared import EGout2, EGout2_2, EGout2_3
thetaarray_sym = np.linspace(-2*10, 2*10, rows*cols)*np.pi/180
m = np.arange(0, int(20), 1)
# EGout2_theta_k=[]
# for i in range(len(thetaarray_sym)):
    # print("##########", i)
N = 100000
xval =2e-3
x = np.linspace(-xval, xval, N)
karray4 = np.linspace(2*np.pi/3e8*(freqi)*1e9, 2*np.pi/3e8*freqf*1e9, N)
EGout2_theta_k = []
for i in range(len(thetaarray_sym)):
    EGout2_theta_ki = EGout2_2(0,0,0,
                L = L, 
                d = d,
                n_g = n_g,
                k = karray4, 
                thetai =thetaarray_sym[i], 
                r = r, 
                m = m, 
                w0 = 8e-4
                )
    EGout2_theta_k.append(EGout2_theta_ki)
    # EGout2_theta_k.append(EGout2i)

#%%
fig, axs = plt.subplots(rows, cols, sharey=False, sharex = False, figsize=(20, 16))

counter = 0
for j in range(cols):
    
    for i in range(rows):
        
        axs[i, j].plot(3e8*karray4*1e-9/(2*np.pi), EGout2_theta_k[counter], color = "r", zorder = len(freq))
        
        axs[rows-1, j].set_xlabel("Frequency (GHz)")
        axs[i, 0].set_ylabel("Amplitude")
        axs[i,j].grid(zorder = 0, alpha = 0.3)
        axs[i,j].set_title(r"$\theta_i$={}°".format(round(thetaarray_sym[counter]*180/np.pi, 5)))
        
        fig.tight_layout(h_pad = 2)
        
        counter += 1
        
# fig.savefig("Gaussian_fixL_fixRho0_VaryTheta,fixWaist1e-3.jpg", dpi = 400)
# fig.savefig("Gaussian_fixL_fixRho0_VaryTheta,fixWaist8e-4.jpg", dpi = 400)

#%%
maxima, _ = find_peaks(EGout2_theta_k[-1], height=[0.00390618340, 0.003906183550])
EGout2_theta_k_maxima = EGout2_theta_k[-1][maxima]
F_maxima = 3e8*karray4[maxima]*1e-9/(2*np.pi)

minima, _ = find_peaks(1/EGout2_theta_k[-1], height=[1/0.003906183425, 1/0.00390618330])
EGout2_theta_k_minima = EGout2_theta_k[-1][minima]
F_minima = 3e8*karray4[minima]*1e-9/(2*np.pi)

m0_maxima = ((EGout2_theta_k_maxima[-1]-EGout2_theta_k_maxima[0]) /
             (F_maxima[-1] - F_maxima[0]))
c0_maxima = EGout2_theta_k_maxima[-1] - m0_maxima*F_maxima[-1]
p0_maxima = [m0_maxima, c0_maxima]
fit_maxima = curve_fit(line, F_maxima, EGout2_theta_k_maxima, p0 = p0_maxima)
p0_maxima_err = np.sqrt(np.diag(fit_maxima[1]))
Farray_maxima = np.linspace(F_maxima[0], F_maxima[-1], 1000)
EGout2_theta_karray_maxima = line(Farray_maxima, *fit_maxima[0])

m0_minima = ((EGout2_theta_k_minima[-1]-EGout2_theta_k_minima[0]) /
             (F_minima[-1] - F_minima[0]))
c0_minima = EGout2_theta_k_minima[-1] - m0_minima*F_minima[-1]
p0_minima = [m0_minima, c0_minima]
fit_minima = curve_fit(line, F_minima, EGout2_theta_k_minima, p0 = p0_minima)
p0_minima_err = np.sqrt(np.diag(fit_minima[1]))
Farray_minima = np.linspace(F_minima[0], F_minima[-1], 1000)
EGout2_theta_karray_minima = line(Farray_minima, *fit_minima[0])

plt.figure()
plt.plot(3e8*karray4*1e-9/(2*np.pi), EGout2_theta_k[-1])
plt.scatter(F_maxima, EGout2_theta_k_maxima, color = "b")
plt.scatter(F_minima, EGout2_theta_k_minima, color = "b")

plt.plot(Farray_maxima, EGout2_theta_karray_maxima, color = "k")
plt.plot(Farray_minima, EGout2_theta_karray_minima, color = "k")
#%%

# =============================================================================
# FINESSE
# =============================================================================

thetaarray_sym2 = np.arange(-10, 10+0.5, 0.5)*np.pi/180
m = np.arange(0, int(20), 1)
# EGout2_theta_k=[]
# for i in range(len(thetaarray_sym)):
    # print("##########", i)
N = 100000
xval =2e-3
x = np.linspace(-xval, xval, N)
freqi2 = freqi
freqf2 = freqi+300
karray5 = np.linspace(2*np.pi/3e8*freqi2*1e9, 2*np.pi/3e8*freqf2*1e9, N)
EGout2_theta_k = []
for i in range(len(thetaarray_sym2)):
    EGout2_theta_ki = EGout2_2(0,0,0,
                L = L, 
                d = d,
                n_g = n_g,
                k = karray5, 
                thetai =thetaarray_sym2[i], 
                r = r, 
                m = m, 
                w0 = 8e-4
                )
    EGout2_theta_k.append(EGout2_theta_ki)
    # EGout2_theta_k.append(EGout2i)
#%%

diffF = []
EGout2_theta_k_intercepts = []
FWHM = []
Finesse = []

for i in range(len(thetaarray_sym2)):
    finessei = finesse(3e8*karray5*1e-9/(2*np.pi), EGout2_theta_k[i], 8, True)
    Finesse.append(finessei)
    
#%%
plt.figure()

plt.scatter(thetaarray_sym2*180/np.pi, Finesse, color = "k", marker = ".", 
            zorder = len(thetaarray_sym2))

plt.grid(zorder=0)

plt.xlabel(r"$\theta_i$°", fontsize = 14)
plt.ylabel("$\mathcal{F}$", fontsize = 14)

plt.xticks(np.arange(-10, 12, 2))
plt.yticks(np.arange(0, 13, 1))

plt.xlim([-10, 10])
plt.ylim([0, 12])

plt.tight_layout()

# plt.savefig("Finesse_vs_tilt_GB.jpg", dpi = 400)

#%%

# =============================================================================
# STUDYING THE EFFECT OF EACH TERM IN IRRADIANCE
# =============================================================================

from Eout_squared import w0w, rho, R, psi, w0wrho
thetaarray_sym = np.linspace(-2*7, 2*8, rows*cols)*np.pi/180
m = np.arange(1, int(500), 1)
# EGout2_theta_k=[]
# for i in range(len(thetaarray_sym)):
    # print("##########", i)
N = 10000
karray4 = np.linspace(2*np.pi/3e8*(freqi)*1e9, 2*np.pi/3e8*freqf*1e9, N)
xval =0.3
x = np.linspace(-xval, xval, N)
y = np.linspace(-xval, xval , N)
# rho = np.linspace(-np.sqrt(x**2+y**2), np.sqrt(x**2+y**2), N)
angle = 13
L=5e-3
w0w_k = w0w(
            x=0,
            y=0,
            z=0,
            L = L, 
            d = d,
            n_g = n_g,
            k = karray4, 
            thetai =angle*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 1e-3
            )

rho_k = rho(x=0,
            y=0,
            z=0,
            L = L, 
            d = d,
            n_g = n_g,
            k = karray4, 
            thetai =angle*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 1e-3)
R_k = R(x=0,
            y=0,
            z=0,
            L = L, 
            d = d,
            n_g = n_g,
            k = karray4, 
            thetai =angle*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 1e-3)

psi_k = psi(x=0,
            y=0,
            z=0,
            L = L, 
            d = d,
            n_g = n_g,
            k = karray4, 
            thetai =angle*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 1e-3)

w0wrho_k = w0wrho(x=0,
            y=0,
            z=0,
            L = L, 
            d = d,
            n_g = n_g,
            k = karray4, 
            thetai =angle*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 1e-3)


#%%  
plt.figure()
plt.plot(np.linspace(freqi, freqf, N), w0w_k)

plt.figure()
plt.plot(np.linspace(freqi, freqf, N), rho_k)

plt.figure()
plt.plot(np.linspace(freqi, freqf, N), R_k)

plt.figure()
plt.plot(np.linspace(freqi, freqf, N), psi_k)
plt.figure()
plt.plot(np.linspace(freqi, freqf, N), w0wrho_k)


#%%
# =============================================================================
# FIX L, FIX RHO = 0, FIX THETA = 0, FIX FREQ = 192300 GHz, VARY WAIST
# =============================================================================

#%%
from Eout_squared import EGout2_2
m = np.arange(0, 69, 1)
w0array = np.linspace(1e-20, 1e-3, int(1e4))
Rarray = np.arange(0.1, 1.0, 0.05)
rarray = np.sqrt(Rarray)
EGout2_w0_r=[]
for i in range(len(rarray)):
    print("##########", i)   
    EGout2i = EGout2_2(0,0,0, 
                L = L, 
                k = 2*np.pi*192300*1e9/3e8, 
                thetai = 0, 
                d = d,
                n_g = n_g,
                r = rarray[i], 
                m = m, 
                w0 = w0array, 
                )
    EGout2_w0_r.append(EGout2i)
    
#%%
plt.figure()
for i in range(len(rarray)):
    plt.plot(w0array, EGout2_w0_r[i], color = "blue", zorder = len(w0array))
    
    plt.xlabel(r"$w_0$")
    plt.ylabel(r"$|E_{out}|^{2}$")
    plt.annotate("{}".format(np.round(Rarray[i], 2)), (w0array[-1]+0.00001, np.max(EGout2_w0_r[i])))
    plt.xlim((w0array[0], w0array[-1]+0.00008))
    plt.ylim((0, 1))
    plt.grid(zorder = 0)
    plt.tight_layout()
plt.annotate("R =", (w0array[-1]+0.00001, np.max(EGout2_w0_r[0])+0.06))
plt.tight_layout()
# plt.savefig("Gaussian_fixLL_fixRho0_fixTheta0,varyWaist, varyr", dpi = 400)

#%%
# =============================================================================
# FIX L, VARY RHO, FIX THETA = 0, FIX FREQ = 192300 GHz, FIX WAIST Many 
# =============================================================================

#%%
from Eout_squared import EGout2_2, psi
m = np.arange(0, int(69), 1)
xval = 4e-3
N = int(1e5)
xarray = np.linspace(-xval, xval, N)
yarray = np.linspace(-xval, xval, N)
rhoarray = np.sqrt(xarray**2 + yarray**2)
rhoarray = np.linspace(-rhoarray[-1], rhoarray[-1], N)
w0array2 = np.linspace(20e-6, 140e-6, rows*cols)

EGout2_rho = []
for i in range(len(w0array2)):
    EGout2_rhoi = EGout2_2(x = xarray, y = yarray, z = 0,
                L = L,
                d = d,
                n_g = n_g,
                k = 2*np.pi*192300*1e9/3e8, 
                thetai = 0*np.pi/180, 
                r = r, 
                m = m, 
                w0 = w0array2[i], 
                )
    EGout2_rho.append(EGout2_rhoi)
    
    
#%%

fig, axs = plt.subplots(rows, cols, sharey=False, sharex = False, figsize=(20, 16))

counter = 0
for j in range(cols):
    
    for i in range(rows):
        
        axs[i, j].plot(rhoarray, EGout2_rho[counter], color = "b", zorder = len(rhoarray), lw = 0.5)
        
        axs[rows-1, j].set_xlabel(r"$\rho$")
        axs[i, 0].set_ylabel(r"$|E_{out}|^{2}$")
        axs[i,j].grid(zorder = 0, alpha = 0.3)
        axs[i,j].set_title(r"$\theta_i$=0, $w_0$={}°".format(round(w0array2[counter], 9)))
        
        fig.tight_layout(h_pad = 2)
        
        counter += 1
        
# plt.savefig("Gaussian_fixLL_varyRho_fixTheta0,fixWaistmany, varyr", dpi = 400)
#%%
xval = 11e-2
N = int(1e5)
xarray = np.linspace(-xval, xval, N)
yarray = np.linspace(-xval, xval, N)
rhoarray = np.sqrt(xarray**2 + yarray**2)
rhoarray = np.linspace(-rhoarray[-1], rhoarray[-1], N)
w0array2 = np.linspace(0.5e-6, 5e-6, rows*cols)

EGout2_rho2 = EGout2_2(x = xarray, y = yarray, z = 0,
            L = L, 
            d = d,
            n_g = n_g,
            k = 2*np.pi*192300*1e9/3e8, 
            thetai = 0*np.pi/180, 
            r = r, 
            m = m, 
            w0 = 8e-2, 
            )
#%%
plt.figure(figsize = (10,6))
# for i in range(1):
plt.plot(rhoarray, EGout2_rho2, color = "blue", zorder = len(w0array))
plt.plot(rhoarray, np.zeros(len(rhoarray))+np.max(EGout2_rho2)/np.exp(2), color = "k", 
         label = r"$|E_{out}|^{2}_{max}/e^2$")
plt.xlabel(r"$\rho$")
plt.ylabel(r"$|E_{out}|^{2}$")
plt.xlim((rhoarray[0], rhoarray[-1]))
plt.xticks(np.arange(-0.15, 0.15+0.025, 0.025))
plt.grid(zorder = 0)
plt.title(r"$w_0$ = 1e-6 m, $\theta^{i}$=0°")
plt.legend()
plt.tight_layout()
# plt.annotate("R =", (w0array[-1]+0.00001, np.max(EGout2_w0_r[0])+0.06))
# plt.tight_layout()

# plt.savefig("Gaussian_fixLL_varyRho_fixTheta0,fixWaistLarge, varyr", dpi = 400)
#%%

from Eout_squared import EGout2_3, reshape_list
m = np.arange(0, int(69), 1)
EGout2_test = EGout2_3(0,0,0,
                L = L, 
                d = d,
                n_g = n_g,
                k = karray4, 
                thetai = 0,
                r = r, 
                m = m, 
                w0 = 1e-3
                )
#%%
plt.figure()
plt.plot(karray4*3e8/(2*np.pi), EGout2_test)
#%%
# =============================================================================
# FIX L, VARY X, FIX THETA = 0, FIX LAMDA = 632.8 NM, FIX WAIST = 1E-3
# =============================================================================

#%%
m = np.arange(0, int(10), 1)


# thetaarray4 = np.zeros(int(1e5))+angle
# xarray = np.linspace(0, 2*69*L*np.sin(angle), len(thetaarray4))
# yarray = np.zeros(len(thetaarray4))

xarray = np.linspace(-0.02, 0.02, int(1e5))
angle = 0*np.pi/180
freq3 = np.linspace(freqi, freqf, 10)*1e9
karray3 = 2*np.pi*freq3/3e8
EGout2_x=[]
EGout2max_x=[]

for i in range(len(thetaarray_sym)):
    EGout2_xi = EGout2_2(x = xarray, y = 0 , z = 0,
                L = L, 
                d = d,
                n_g = n_g,
                k = karray[0], 
                thetai =  thetaarray_sym[i], 
                r = r, 
                m = m, 
                w0 = 8e-4, 
                )
    EGout2_x.append(EGout2_xi)
    EGout2max_x.append(np.max(EGout2_xi))
    

#%%
fig, axs = plt.subplots(rows, cols, sharey=False, sharex = False, figsize=(20, 16))

counter = 0
for j in range(cols):
    
    for i in range(rows):
        
        axs[i, j].plot(xarray, EGout2_x[counter], color = "r", zorder = len(freq))
        
        axs[rows-1, j].set_xlabel("Frequency (GHz)")
        axs[i, 0].set_ylabel("Amplitude")
        axs[i,j].grid(zorder = 0, alpha = 0.3)
        axs[i,j].set_title(r"$\theta_i$={}°".format(round(thetaarray_sym[counter]*180/np.pi, 5)))
        
        fig.tight_layout(h_pad = 2)
        
        counter += 1
        
# fig.savefig("Gaussian_varyx_fixL_fixF_fixRho0_VaryTheta_fixWaist8e-4.jpg", dpi = 400)



