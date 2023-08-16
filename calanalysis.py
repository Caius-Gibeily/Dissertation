#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 20:44:35 2023

@author: caiusgibeily
"""


#### Aethon and Ca2+ analysis ####
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
os.chdir("/home/caiusgibeily/Documents/Oxdrive/Project_2/Dissertation/Figure 10 - Aethon/Data/")

def plot_pulses(dat, sfactor):
    plt.figure(figsize=(20,12))
    x = np.arange(0,(len(dat[10:])*sfactor)/100,sfactor/100)
    plt.plot(x,dat["Mean"].iloc[10:],linewidth=2)
    plt.xlabel("Time (s)",fontsize=40)
    plt.ylabel("Intensity (a.u.)",fontsize=40)
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    
    
two_hz = pd.read_csv("2hz.csv",sep=",")
five_hz = pd.read_csv("5hz - 3.csv",sep=",")
ten_hz = pd.read_csv("10hz- 1.csv",sep=",")
five_hz_other = pd.read_csv("5hz - other 3.csv",sep=",")
##
os.chdir("/home/caiusgibeily/Documents/Oxdrive/Project_2/Dissertation/Figure 10 - Aethon/Data/low_int")

two_hza = pd.read_csv("2hz-1.csv",sep=",")
five_hza = pd.read_csv("5hz-1.csv",sep=",")
ten_hza = pd.read_csv("10hz-1.csv",sep=",")
five_hz_other = pd.read_csv("10hz-1-other.csv",sep=",")
##

plot_pulses(two_hz,5)
plot_pulses(five_hz,3)
plot_pulses(ten_hz,1)
plot_pulses(five_hz_other.iloc[60:],3)
plt.xlim([0,10])
###
plot_pulses(two_hza,5)
plot_pulses(five_hza,3)
plot_pulses(ten_hza,1)
plot_pulses(five_hz_other.iloc[60:],3)
plt.xlim([0,10])



close_ten = plot_pulses(ten_hz.iloc[20:200],1)

##
os.chdir("/home/caiusgibeily/Documents/Oxdrive/Project_2/Dissertation/Figure 12 - Ca analysis")

gcamp = pd.read_csv("fluo4.csv",sep=",")
plt.figure(figsize=(18,15))
for i in range(len(gcamp.iloc[0,:])//10):
    #if gcamp["Mean" + str(i+1)].iloc[-1] - gcamp["Mean" + str(i+1)].iloc[0] > 1:
    plt.plot((gcamp["Mean" + str(i+1)]-gcamp["Mean" + str(i+1)].iloc[0])/gcamp["Mean" + str(i+1)].iloc[0])
    
plt.ylim([-0.2,0.2])
plt.ylabel("FLuorescence (ΔF/F)",fontsize=40)
plt.xlabel("Time (frames)",fontsize=40)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.plot(gcamp["Mean" + str(i+1)].iloc[0:1000])

## Transfection dat
## mean by row

col = "SH - Cell Length [µm] - Mean per Well"
trans = pd.read_csv("trans_data.csv",sep=",")
means = trans.groupby(["Row"])[col].mean().reset_index().pivot_table(index="Row", values=col,fill_value=-0.1)
stds = trans.groupby(["Row"])[col].std().reset_index().pivot_table(index="Row", values=col,fill_value=-0.1)/np.sqrt(2)

plt.figure(figsize=(18,15))
plt.scatter(trans["Row"],trans[col],s=100)
plt.errorbar([0.1,0.25,0.5,0.75,1],means[col],stds[col],elinewidth=8,color="black",capthick=5,capsize=10,linewidth=6)
plt.xlabel("Stock viral volume (μl)",fontsize=40)
plt.ylabel("Cell Length (μm)",fontsize=40)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.xlim([1.1,0.08])

from sklearn.linear_model import LinearRegression, LogisticRegression
from scipy.special import expit
X = np.array(trans[col])[:, np.newaxis]
y = np.array(trans["Row"]).astype(float)

