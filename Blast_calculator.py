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
blast_strenght = 10

def sorting(arr, var):
    pos = bisect.bisect_left(arr, var)
    if pos == 0:
        return arr[0]
    if pos == len(arr):
        return arr[-1]
    before = arr[pos-1]
    after = arr[pos]
    return [after, before]

def energy(mass=None):
    if mass==None:
        mass=1
    volume = mass/(0.08988*29.5)
    return 3500000*volume

def scaled_dist(dist=None):
    if dist==None:
        n = 100/((P_a/energy())**(1/3))
        distance = np.arange(1, n+1, 1, dtype=float)
    scaled_distance_coeff = ((P_a/energy())**(1/3))
    scaled_distance = np.multiply(dist, scaled_distance_coeff)
    return scaled_distance

def overpressure(dist=None):
    P_file = json.load(
        open(path + "\\overpressure\\{}.json".format(blast_strenght), 'r'))
    P_df = pd.DataFrame(P_file)
    P_arr = P_df.to_numpy()
    P_arr=P_arr.astype(np.float32)

    if dist == None:
        k = np.interp(scaled_dist(), P_arr[:, 0], P_arr[:, 1])
        overpressure = np.multiply(k, 100)
    else:
        s = sorting(P_arr[:,0], dist)
        ar0 = P_arr[:, 0]
        ar1 = P_arr[:, 1]
        v=[]
        for i in s:
            ind=np.where(ar0 == i)
            v.append(ar1[ind])
        k = v[0] + (scaled_dist(dist)- s[0])* (v[1]-v[0])/(s[1]-s[0])
        overpressure = np.multiply(k, 100)
    return overpressure


def duration(dist=None):
    T_file = json.load(
        open(path + "\\duration\\{}.json".format(blast_strenght), 'r'))
    T_df = pd.DataFrame(T_file)
    T_arr = T_df.to_numpy()
    T_arr = T_arr.astype(np.float32)

    calc1 = ((energy()/P_a) ** (1/3)) * 1/a
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

def plotting():
    fig, axs = plt.subplots(3, sharex = True, sharey = False)
    axs[0].plot(scaled_dist(), overpressure(), '-')
    axs[1].plot(scaled_dist(), duration(), '-')
    axs[2].plot(scaled_dist(), impulse(), '-')
    axs[0].set_title('side-on overpressure, kPa')
    axs[1].set_title('positive phase duration, s')
    axs[2].set_title('impulse, kPa*s')
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    plt.show()
    
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
    pb_file = json.load(open(path + "probit.json", 'r'))
    pb_df = pd.DataFrame(pb_file)
    pb_arr = pb_df.to_numpy()
    pb_arr = pb_arr.astype(np.float32)
    prb = np.interp(probit, pb_arr[:, 0], pb_arr[:, 1])

if __name__ == '__main__':
    m = input("Enter mass of hydrogen: ")
    if m.isnumeric()==True:
        mass = int(m)
    elif m == None:
        mass = 1
    else:
        print("Wrong input")
        sys.exit(0)
    d = input("Enter distance between wall and explosion center: ")
    if d.isnumeric() == True:
        dist = int(d)
        pb1 = probit_for_structure_damage(dist)
        pb2 = window_breakage(dist)
        pb3 = frag_velocity(dist)
        #perc = 
        print(pb1)
        print(pb2)
        print(pb3)
    elif d == None:
        scaled_dist()
        overpressure()
        impulse()
        plotting()
    else:
        print("Wrong input")
        sys.exit(0)
    
'''    parser = argparse.ArgumentParser(description='Process initial values')
    parser.add_argument('-m', '--mass', type=int, metavar='Mass of hydrogen', default=1, help='Mass of hydrogen in stoichiometric hydrogen-air mixture')
    parser.add_argument('-w', '--wall', type=int, metavar='Distance to obstacle', default=1,
                        help='Distance between center of explosion and obstacle')
    parser.add_argument('-p', '--plot', action='store_true', default=False,
                        help='Plot blast parameters') '''
