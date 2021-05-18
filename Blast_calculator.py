import math
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse
import json
import sys
import bisect

atmospheric_pressure_default = 101250
a = 340 #m/s
P_a = 101325 #Pa
path = os.getcwd()

def sorting(arr, var):
    pos = bisect.bisect_left(arr, var)
    if pos == 0:
        return arr[0]
    if pos == len(arr):
        return arr[-1]
    before = arr[pos-1]
    after = arr[pos]
    return [after, before]

def energy(mass):
    volume = mass*11.2*100/29.5
    return 3500000*volume

def scaled_dist(dist=None):
    if dist==None:
        n = 100/((P_a/energy(mass))**(1/3))
        dist = np.arange(1, n+1, 1, dtype=float)
    scaled_distance_coeff = ((P_a/energy(mass))**(1/3))
    scaled_distance = np.multiply(dist, scaled_distance_coeff)
    return scaled_distance

def overpressure(dist=None):
    P_file = json.load(
        open(path + "\\overpressure\\{}.json".format(blast_strength), 'r'))
    P_df = pd.DataFrame(P_file)
    P_arr = P_df.to_numpy()
    P_arr=P_arr.astype(np.float32)

    if dist == None:
        k = np.interp(scaled_dist(), P_arr[:, 0], P_arr[:, 1])
        overpressure = np.multiply(k, P_a)
    else:
        s = sorting(P_arr[:,0], dist)
        ar0 = P_arr[:, 0]
        ar1 = P_arr[:, 1]
        v=[]
        for i in s:
            ind=np.where(ar0 == i)
            v.append(ar1[ind])
        k = v[0] + (scaled_dist(dist)- s[0])* (v[1]-v[0])/(s[1]-s[0])
        overpressure = np.multiply(k, P_a)
    return overpressure

def duration(dist=None):
    T_file = json.load(
        open(path + "\\duration\\{}.json".format(blast_strength), 'r'))
    T_df = pd.DataFrame(T_file)
    T_arr = T_df.to_numpy()
    T_arr = T_arr.astype(np.float32)

    calc1 = ((energy(mass)/P_a) ** (1/3)) * 1/a
    if dist == None:
        l = np.interp(scaled_dist(), T_arr[:, 0], T_arr[:, 1])
        positive_phase_duration = np.multiply(calc1, l)
    else:
        s = sorting(T_arr[:, 0], dist)
        ar0 = T_arr[:, 0]
        ar1 = T_arr[:, 1]
        u=[]
        for i in s:
            ind=np.where(ar0 == i)
            u.append(ar1[ind])
        l = u[0] + (scaled_dist(dist) - s[0]) * (u[1]-u[0])/(s[1]-s[0])
        positive_phase_duration = np.multiply(calc1, l)
    return positive_phase_duration

def impulse(dist=None):
    impulse = np.multiply(0.5*overpressure(dist), duration(dist))
    return impulse

def plotting(x, y, ax=None, plt_kwargs={}):
    if ax is None:
        ax = plt.gca()
    ax.plot(x,y,**plt_kwargs)
    return(ax)
    
def probit_for_structure_damage(dist):
    sd = scaled_dist(dist)
    ovp = overpressure(sd)
    imp = impulse(sd)
    if ovp <= 4600 and imp <= 110:
        V = ((4600/ovp)**3.9)+((110/imp)**5)
        return 5 - 0.26*np.log(V)
    elif ovp >= 40000 and imp >= 460:
        V = ((40000/ovp)**7.4)+((460/imp)**11.3)
        return 5 - 0.22*np.log(V)
    else:
        V = ((17500/ovp)**8.4)+((290/imp)**9.2)
        return 5 - 0.26*np.log(V)

def window_breakage(dist):
    sd = scaled_dist(dist)
    ovp = overpressure(sd)
    return -11.97+2.12* np.log(ovp)

def frag_velocity(dist):
    sd = scaled_dist(dist)
    ovp = overpressure(sd)
    P_r = 2*ovp + ((2.4*ovp**2)/(0.4*ovp)+2*1.4*P_a)
    
    return P_r

def percentage(probit):
    pb_file = json.load(open(path + "\\probit.json", 'r'))
    pb_df = pd.DataFrame(pb_file)
    pb_arr = pb_df.to_numpy()
    pb_arr = pb_arr.astype(np.float32)
    prb = np.interp(probit, pb_arr[:, 0], pb_arr[:, 1])
    return prb

if __name__ == '__main__':
    m = input("Enter mass of hydrogen: ")
    if m.isnumeric()==True:
        mass = int(m)
    else:
        print("Wrong input. Mass will be set as 1 kg by default")
        mass = 1
    blast_strength = input("Enter blast strength (from 1 to 10): ")
    if blast_strength.isnumeric() == False and blast_strength > 10:
        print("Wrong input. Blast strength will be set to 10 by default")
        blast_strength = 10        
    else:
        pass
    d = input("Enter distance between wall and explosion center: ")
    if d.isnumeric() == True:
        dist = int(d)
        pb1 = round(percentage(float(probit_for_structure_damage(dist))), 3)
        if pb1<=99.9:
            print(f'The probability of damage to building itself is {pb1}%')
        elif pb1>99.9:
            print(f'The probability of damage to building itself is higher than 99.9%')
        pb2 = round(percentage(float(window_breakage(dist))),3)
        if pb1 <= 99.9:
            print(f'The probability of windows breakage is {pb2}%')
        elif pb1 > 99.9:
            print(f'The probability of windows breakage is higher than 99.9%')
        pb3 = frag_velocity(dist)

        print(pb3)
    else:
        print("Wrong input")
        pass
    pl=input("Are plots of overpressure and impulse necessary? y/n ")
    if pl.lower == 'y':
        x = scaled_dist()/((P_a/energy(mass))**(1/3))
        plotting(x, overpressure()/1000, title='side-on overpressure, kPa')
        plt.plot(x, duration(), title='positive phase duration, s')
        plt.plot(x, impulse()/1000, title='impulse, kPa*s')
        plt.grid()
        plt.show()
    elif pl.lower == 'n':
        sys.exit(0)
