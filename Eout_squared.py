# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 23:23:38 2022

@author: pgara
"""

#editing

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt



def reshape_list(oglist, nest_len):
    
    new_list = [oglist[i:i+nest_len] for i in range(0, len(oglist), nest_len)]
    new_list = [l.tolist() for l in new_list]

    return new_list

def Eout2_norm(L, k, r1, r2, A=1):
    
    return A**2*(1-r1**2)*(1-r2**2)/(r1**2*r2**2 - 2*r1*r2*np.cos(2*k*L) +1)


def Eout2_tilt(L, k, thetai, r, A=1):
    
    L1 = L/np.cos(thetai)
    L2 = L1*abs(np.cos(2*thetai))
    
    t = np.sqrt(1-r**2)
    
    return A**2*t**4/(r**4 - 2*r**2*np.cos(k*(L1+L2)) +1)

def FSR_tilt(L, thetai):
    c=3e8
    L1 = L/np.cos(thetai)
    L2 = L1*abs(np.cos(2*thetai))
    return c/(L1+L2)

def EGout2(x, y, L, d, n_g, k, thetai, r, m, w0, A=1):
    
    n_vac = 1
    t = np.sqrt(1-r**2)
    
    L1 = L/np.cos(thetai)
    L2 = L1*abs(np.cos(2*thetai))
    zm = (m+1)*L1 + m*L2
    
    xm = x-2*m*L*np.sin(thetai)
    rhom = np.sqrt(xm**2 + y**2)
    lamda = 2*np.pi/k 
    zR = np.pi*w0**2/lamda
    w = w0*np.sqrt(1 + (zm/zR)**2)
    psi = np.arctan(zm/zR)
    R = zm + zR**2/zm
    
    EGm = (t**4*r**(2*m) * 
            A * 
            w0/w * 
            np.exp(-rhom**2/w**2) * 
            np.exp(-1j*k*zm) * 
            np.exp(-1j*k*rhom**2/(2*R)) * 
            np.exp(1j*psi))
    EGmstar = np.conj(EGm)
    
    Eout = sum(EGm)
    Eoutstar = sum(EGmstar)
    
    return np.real(Eout*Eoutstar)

def EGout2(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    
    t = np.sqrt(1-r**2)
    n_vac = 1
    thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
    L1 = L/np.cos(thetai)
    L2 = L1*abs(np.cos(2*thetai))
    D = d/np.cos(thetar)*np.cos(thetai-thetar)
    zm = z + (m+1)*L1 + m*L2 + 2*D
    
    xm = x-2*m*L*np.sin(thetai)
    rhom = np.sqrt(xm**2 + y**2)
    lamda = 2*np.pi/k 
    zR = np.pi*w0**2/lamda
    w = w0*np.sqrt(1 + (zm/zR)**2)
    psi = np.arctan(zm/zR)
    R = zm + zR**2/zm
    
    EGm = (t**4*r**(2*m) * 
            A * 
            w0/w * 
            np.exp(-rhom**2/w**2) * 
            np.exp(-1j*k*zm) * 
            np.exp(-1j*k*rhom**2/(2*R)) * 
            np.exp(1j*psi))
    EGmstar = np.conj(EGm)
    
    Eout = sum(EGm)
    Eoutstar = sum(EGmstar)
    
    return np.real(Eout*Eoutstar)


def PGout(rho, L, k, thetai, r, m, w0, A=1):
    
    return EGout2(rho, L, k, thetai, r, m, w0, A)*2*np.pi*rho

def EGout2_2(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    Eoutstar = 0
    
    for M in m:
        print("m", M)
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm = (t**4*r**(2*M) * 
                A * 
                w0/w * 
                np.exp(-rhom**2/w**2) *
                np.exp(-1j*k*zm) *
                np.exp(-1j*k*rhom**2/(2*R))*
                np.exp(1j*psi))
                
        EGmstar = np.conj(EGm)
        Eout += EGm
        
        # print(w0/w)
        # print(np.exp(-rhom**2/w**2))
        # print(w0/w*np.exp(-rhom**2/w**2))
        
        Eoutstar += EGmstar
    
    return np.real(Eout*Eoutstar)

def PGout_2(rho, x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    
    Eout = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        print(xm)
        # rhom = np.sqrt(xm**2 + y**2)        
        # rho = np.sqrt(x**2+y**2)
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm = (t**4*r**(2*M) * 
                A * 
                w0/w * 
                np.exp(-rho**2/w**2) *
                np.exp(-1j*k*zm))
                # np.exp(-1j*k*rhom**2/(2*R))*
                # np.exp(1j*psi))
                
        # EGmstar = np.conj(EGm)
        Eout += EGm
        # Eoutstar += EGmstar
    
    return Eout*np.conj(Eout)

def EGout2_3(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    
    
    m = np.repeat(m, len(k))
    m = np.array(reshape_list(m, len(k)))
    t = np.sqrt(1-r**2)
    n_vac = 1
    thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
    L1 = L/np.cos(thetai)
    L2 = L1*abs(np.cos(2*thetai))
    D = d/np.cos(thetar)*np.cos(thetai-thetar)
    zm = z + (m+1)*L1 + m*L2 + 2*D
    
    xm = x-2*m*L*np.sin(thetai)
    rhom = np.sqrt(xm**2 + y**2)
    lamda = 2*np.pi/k 
    zR = np.pi*w0**2/lamda
    w = w0*np.sqrt(1 + (zm/zR)**2)
    psi = np.arctan(zm/zR)
    R = zm + zR**2/zm
    
    EGm = (t**4*r**(2*m) * 
            A * 
            w0/w * 
            np.exp(-rhom**2/w**2) * 
            np.exp(-1j*k*zm) * 
            np.exp(-1j*k*rhom**2/(2*R)) * 
            np.exp(1j*psi))
    EGmstar = np.conj(EGm)
    
    Eout = sum(EGm)
    Eoutstar = sum(EGmstar)
    
    return np.real(Eout*Eoutstar)


def w0w(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    # Eoutstar = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x-2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm =  (t**4*r**(2*M) * 
                w0/w*
                np.exp(-1j*k*zm))
                # np.exp(-1j*k*rhom**2/(2*R))*
                # np.exp(1j*psi))
                
        # EGmstar = np.conj(EGm)
        print(zm)
        print(zR)
        Eout += EGm
        
        
        # Eoutstar += EGmstar
    
    return Eout*np.conj(Eout)  

def rho(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    # Eoutstar = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm =  (t**4*r**(2*M) * 
        
                np.exp(-1j*k*zm)*
                np.exp(-rhom**2/w**2))
                # np.exp(1j*psi))
                
        # EGmstar = np.conj(EGm)
        Eout += EGm
        
        
        # Eoutstar += EGmstar
    
    return Eout*np.conj(Eout) 

def R(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    # Eoutstar = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm =  (t**4*r**(2*M) * 
                np.exp(-1j*k*rhom**2/(2*R))*
                np.exp(-1j*k*zm))
                # np.exp(-rhom**2/w**2))
                # np.exp(1j*psi))
                
        # EGmstar = np.conj(EGm)
        Eout += EGm
        
        
        # Eoutstar += EGmstar
    
    return Eout*np.conj(Eout) 

def psi(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    Eoutstar = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm =  (t**4*r**(2*M) * 
                # np.exp(-1j*k*rhom**2/(2*R)))
                np.exp(-1j*k*zm)*
                # np.exp(-rhom**2/w**2))
                np.exp(1j*psi))
                
        EGmstar = np.conj(EGm)
        Eout += EGm
        
        
        Eoutstar += EGmstar
    
    return Eout*Eoutstar

def w0wrho(x, y, z, L, d, n_g, k, thetai, r, m, w0, A=1):
    """
    Same as the previous one but iterates over each m element. Useful when k has a
    different shape than m.
    """
    Eout = 0
    # Eoutstar = 0
    
    for M in m:
    
        t = np.sqrt(1-r**2)
        n_vac = 1
        thetar = np.arcsin(n_vac/n_g*np.sin(thetai))
        L1 = L/np.cos(thetai)
        L2 = L1*np.cos(2*thetai)
        D = d/np.cos(thetar)*np.cos(thetai-thetar)
        zm = z + (M+1)*L1 + M*L2 +2*D
        # print(zm)

        xm = x+2*M*L*np.sin(thetai)
        rhom = np.sqrt(xm**2 + y**2)        
        
        lamda = 2*np.pi/k 
        # lamda = np.sort(lamda)
        # print(lamda)
        zR = np.pi*w0**2/lamda
        w = w0*np.sqrt(1+(zm/zR)**2)
        # print(w)
        psi = np.arctan(zm/zR)
        R = zm + zR**2/zm
        
        EGm =  (t**4*r**(2*M) * 
                w0/w *
                np.exp(-1j*k*rhom**2/(2*R))*

                np.exp(-1j*k*zm)*
                np.exp(-rhom**2/w**2)*
                np.exp(1j*psi))
                
        # EGmstar = np.conj(EGm)
        Eout += EGm
        
        
        # Eoutstar += EGmstar
    
    return Eout*np.conj(Eout) 

def line(F, m, c):
    
    return F*m +c

#_____________________________________________________________________________

def finesse(F, V, acc, display = False):
    """
    acc: precision to which the frequencies corresponding to the max amplitude/2
    are obtained
    """
    hrange = np.max(V)-np.min(V)
    
    # Find maxima
    hmin_maxima = np.max(V)-hrange/2
    hmax_maxima = np.max(V)+hrange/2
    maxima, _ = find_peaks(V, height=[hmin_maxima, hmax_maxima])
    V_maxima = V[maxima]
    F_maxima = F[maxima]
    
    # Find minima
    hmin_minima = np.min(V)-hrange/5
    hmin_minima = 1e-6
    print("hmin_minima", hmin_minima)
    hmax_minima = np.min(V)+hrange/2
    print("hmax_minima", hmax_minima)
    minima, _ = find_peaks(1/V, height=[1/hmax_minima, 1/hmin_minima])
    V_minima = V[minima]
    F_minima = F[minima]
    
    # Fit maxima
    m0_maxima = ((V_maxima[-1]-V_maxima[0]) /
                 (F_maxima[-1] - F_maxima[0]))
    c0_maxima = V_maxima[-1] - m0_maxima*F_maxima[-1]
    p0_maxima = [m0_maxima, c0_maxima]
    fit_maxima = curve_fit(line, F_maxima, V_maxima, p0 = p0_maxima)
    p0_maxima_err = np.sqrt(np.diag(fit_maxima[1]))
    Farray_maxima = np.linspace(F_maxima[0], F_maxima[-1], 1000)
    Varray_maxima = line(Farray_maxima, *fit_maxima[0])
    
    # Fit minima
    m0_minima = ((V_minima[-1]-V_minima[0]) /
                 (F_minima[-1] - F_minima[0]))
    c0_minima = V_minima[-1] - m0_minima*F_minima[-1]
    p0_minima = [m0_minima, c0_minima]
    fit_minima = curve_fit(line, F_minima, V_minima, p0 = p0_minima)
    p0_minima_err = np.sqrt(np.diag(fit_minima[1]))
    Farray_minima = np.linspace(F_minima[0], F_minima[-1], 1000)
    Varray_minima = line(Farray_minima, *fit_minima[0])
    
    # Intercepts
    m_avg = (fit_maxima[0][0] + fit_minima[0][0])/2
    c_avg = (fit_maxima[0][1] + fit_minima[0][1])/2
    midline = line(F, m_avg, c_avg)
    print(midline)
    
    intercept_indxs = np.argwhere(np.diff(np.sign(midline - V))).flatten()
    
    Fintercepts = F[intercept_indxs]
    Vintercepts = V[intercept_indxs]
    if F[0] < F[maxima[0]] and V[0] > midline[0]:
        Fintercepts = np.delete(Fintercepts, 0)
        Vintercepts = np.delete(Vintercepts, 0)

    if F[maxima[-1]] > F[-1] and V[-1] < midline[-1]:
        Fintercepts = np.delete(Fintercepts, -1)
        Vintercepts = np.delete(Vintercepts, -1)
        
    

    Fintercepts2 = np.delete(Fintercepts, 0)
    Fintercepts2 = np.insert(Fintercepts2, len(Fintercepts2), 0)
        
    diffFintercepts = np.delete(Fintercepts2-Fintercepts, -1)
    FWHM = diffFintercepts[0::2]
    FWHM_avg = sum(FWHM)/len(FWHM)
    
    F_maxima2 = np.delete(F_maxima, 0)
    F_maxima2 = np.insert(F_maxima2, len(F_maxima2), 0)
    
    FSR = np.delete(F_maxima2-F_maxima, -1)
    # print(FSR, len(FSR))
    FSR_avg = sum(FSR)/len(FSR)
    
    finesse_avg = FSR_avg/FWHM_avg
    # Plot
    
    if display == True:
        plt.figure()
        plt.plot(F, V, color = "r")
        plt.scatter(F_maxima, V_maxima, color = "b")
        plt.scatter(F_minima, V_minima, color = "b")
        plt.scatter(Fintercepts, Vintercepts, color = "k", marker = ".", zorder = 2000)
    
        plt.plot(Farray_maxima, Varray_maxima, color = "k")
        plt.plot(Farray_minima, Varray_minima, color = "k")
        plt.plot(F, midline, color = "lime")
        
        plt.xlabel("Frequency (GHz)")
        plt.ylabel(r"$|E_{out}|^{2}$")
        plt.title(r"$\theta_i$=18°")
        plt.tight_layout()
        
        
    return finesse_avg
    
    
    
    
    
    
        
        
    
    
    
    
    
