#! / usr / bin /env python3

import math
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse
import json

atmospheric_pressure_default = 101250
a = 340 #m/s
P_a = 101325 #Pa
path = os.getcwd()
blast_strenght = 10
mass = 1 #kg
volume = mass/(0.08988*29.5)
E = 3500000*volume #J

def blast_params():
    P_file = json.load(
        open(path + "\\overpressure\\{}.json".format(blast_strenght), 'r'))
    P_df = pd.DataFrame(P_file)
    P_arr = P_df.to_numpy()
    P_arr=P_arr.astype(np.float32)

    T_file = json.load(
        open(path + "\\duration\\{}.json".format(blast_strenght), 'r'))
    T_df = pd.DataFrame(T_file)
    T_arr = T_df.to_numpy()
    T_arr = T_arr.astype(np.float32)

    n = 100/((P_a/E)**(1/3))
    distance = np.arange(1, n+1, 1, dtype=float)
    scaled_distance_coeff = ((P_a/E)**(1/3))
    scaled_distance = np.multiply(distance, scaled_distance_coeff)
    k = np.interp(scaled_distance, P_arr[:, 0], P_arr[:, 1])
    #add 0.24
    overpressure = np.multiply(k, 100)

    l = np.interp(scaled_distance, T_arr[:, 0], T_arr[:,1])
    calc1 = ((E/P_a) ** (1/3)) * 1/a
    positive_phase_duration = np.multiply(calc1, l)

    impulse = np.fromfunction(
        lambda overpressure, positive_phase_duration: 0.5*overpressure*positive_phase_duration, (1, len(distance)+1), dtype=float)
    impulse = np.multiply(0.5*overpressure, positive_phase_duration)

#print(len(overpressure))
#print(T_arr)
#print(P_arr)
#print(len(positive_phase_duration))
#print(len(impulse))
#print(len(scaled_distance))
#print(o)
#print(T_s)
#print(P_s)
#print(dist)
#print(to_add)
#print(volume)
def plotting_blast():
    fig, axs = plt.subplots(3, sharex = True, sharey = False)
    axs[0].plot(scaled_distance, overpressure, '-')
    axs[1].plot(scaled_distance, positive_phase_duration, '-')
    axs[2].plot(scaled_distance, impulse, '-')
    axs[0].set_title('side-on overpressure, kPa')
    axs[1].set_title('positive phase duration, s')
    axs[2].set_title('impulse, kPa*s')
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    plt.show()

def 

parser = argparse.ArgumentParser(description='Process initial values.')
parser.add_argument('Mass', )
