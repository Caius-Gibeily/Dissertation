#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:10:40 2023

@author: caiusgibeily
"""

############################################
########## Image processing ################
############################################
import numpy as np
import os, glob
from PIL import Image
import natsort
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm
from  matplotlib.colors import LinearSegmentedColormap
MSN2 = "/run/media/caiusgibeily/NEW EVOS DRIVE/cg123/MSN_plate2/MSN_plate2__2023-07-15T16_51_57-Measurement 1/Images"
DA2 = "/run/media/caiusgibeily/NEW EVOS DRIVE/cg123/DA/DA2__2023-07-30T19_31_27-Measurement 2/Images"
imdir = DA2
#savedir = "/home/caiusgibeily/Documents/Oxdrive/Project_2/Im_analysis/"
savedir = "/run/media/caiusgibeily/NEW EVOS DRIVE/cg123/DA/DA-split/"
os.chdir(imdir)

### iterator ##
wells_MSN2 = ["B5","B6","B7","B8","B9",
         "C5","C6","C7","C8","C9",
         "D5","D6","D7","D8","D9",
         "E5","E6","E7","E8","E9",
         "F5","F6","F7","F8","F9",
         "G5","G6","G7","G8","G9"]
wells_DA = ["B2","B3","B4","B5","B6","B7",
         "C2","C3","C4","C5","C6","C7",
         "D2","D3","D4",
         "E2","E3","E4",
         "F2","F3","F4",
         "G2","G3"]


wells = wells_DA

images = glob.glob("*.tiff") 
images = natsort.natsorted(images)
zpos = 4
FOVs = 25

os.chdir(savedir)
save_ims = glob.glob("*.tiff")
save_im = natsort.natsorted(save_ims)

channels = ["blue","red","yellow","green"] #
n_channels = len(channels)
well_elements = FOVs*n_channels*zpos

file_list = np.array([save_ims[i][0:4] for i in range(len(save_ims))]).astype(int)
if len(file_list) != 0:
    iterator = max(file_list)
else:
    iterator = 0
os.chdir(imdir)
for i,well in enumerate(wells):
    i = int(iterator/100)
    well = wells[i]
    im_range = images[i*well_elements:(i+1)*well_elements]
    for j in range(FOVs):
        FOV_range = im_range[j*n_channels*zpos:(j+1)*n_channels*zpos]
        for k, channel in enumerate(channels):
            for l in range(zpos):
                im = Image.open(FOV_range[k+n_channels*l])
                #array = np.asarray(im, dtype=np.int32)
                #stack = np.dstack((stack,array))
            #mean = np.mean(stack[:,:,3:], axis=2)
            #mean = stack[:,:,2]
            # Create a new grayscale image from the grayscale value
            #mean_img = Image.fromarray(mean)
            im.save(savedir + str(iterator).zfill(4) + "_" + well + "_" + str(j).zfill(2) + "_" + channel + ".tiff", "TIFF")
            #mean_img.save(savedir + str(iterator) + ".tiff", "TIFF")
            
            iterator += 1
    
####
# Image colocalisation analyses #

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage as ndi
from skimage import data, filters, measure, segmentation
import pandas as pd
import math

#############################
#### Intensity plot #########
#############################
def jitter(arr, stdev = 0.06):
    return arr + np.random.randn(len(arr)) * stdev

###
imdir = "/home/caiusgibeily/Documents/Oxdrive/Project_2/Im_analysis/"
os.chdir(imdir)

########################### MSNs ############################################## 
intensities = pd.read_csv("intensities_MSN.csv")
antibodies = np.array(intensities.iloc[:,1:4]) ### first four columns contain antibody data
antibodies = np.unique(antibodies); antibodies = antibodies[antibodies!="0"]
## aggregate intensities by receptor
channels = 3
aggregated = {}
plt.figure(figsize=(20,15))
for i, anti in enumerate(antibodies):    
    aggregated[anti] = np.zeros((1,4))
    for j in range(channels):
        rec_dat = intensities.loc[intensities.iloc[:,1+j] == anti]
        if len(rec_dat)!= 0:
            rec_dat = rec_dat.iloc[:,4*j+6:4*j+10]
            aggregated[anti] = np.append(aggregated[anti],rec_dat,0)
            aggregated[anti] = aggregated[anti][1:,:]
    
    normalised = aggregated[anti][:,0]/(aggregated[anti][:,2]-aggregated[anti][:,3])
    x = jitter(np.repeat(i+1,len(aggregated[anti])),0.1)
    plt.errorbar(x, normalised, yerr=aggregated[anti][:,1]/(aggregated[anti][:,2]-aggregated[anti][:,3]), fmt="o",c="black")
    plt.scatter(x, normalised,c="black",s=40)
    #plot_array= np.append(np.full((len(aggregated[anti]),1),i),np.reshape(normalised,(len(normalised),1)),1)
data = [aggregated[anti][:,0]/(aggregated[anti][:,2]-aggregated[anti][:,3]) for i,anti in enumerate(antibodies)]
violin = plt.violinplot(data,showmedians=True)
for pc in violin["bodies"]:
    pc.set_facecolor("red")
    pc.set_edgecolor("black")
    pc.set_alpha(0.5)
for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    vp = violin[partname]
    vp.set_edgecolor("black")
    vp.set_linewidth(3)
plt.xticks(np.arange(1,len(antibodies)+1,1),antibodies,fontsize=40)
plt.yticks(fontsize=40)
plt.xlabel("Protein",fontsize=40)
plt.ylabel("Normalised intensity",fontsize=40)


######################### DANs ###################################
intensities = pd.read_csv("intensities_DAN.csv")

antibodies = np.array(intensities.iloc[:,1:3])
antibodies = np.unique(antibodies); antibodies = antibodies[antibodies!="0"]
## aggregate intensities by receptor
channels = 3
aggregated = {}
plt.figure(figsize=(20,15))
for i, anti in enumerate(antibodies):    
    aggregated[anti] = np.zeros((1,4))
    for j in range(channels):
        rec_dat = intensities.loc[intensities.iloc[:,1+j] == anti]
        if len(rec_dat)!= 0:
            rec_dat = rec_dat.iloc[:,4*j+5:4*j+9]
            aggregated[anti] = np.append(aggregated[anti],rec_dat,0)
            aggregated[anti] = aggregated[anti][1:,:]
    
    normalised = aggregated[anti][:,0]/(aggregated[anti][:,2])
    x = jitter(np.repeat(i+1,len(aggregated[anti])),0.1)
    plt.errorbar(x, normalised, yerr=aggregated[anti][:,1]/aggregated[anti][:,2], fmt="o",c="black")
    plt.scatter(x, normalised,c="black",s=40)
    #plot_array= np.append(np.full((len(aggregated[anti]),1),i),np.reshape(normalised,(len(normalised),1)),1)
data = [aggregated[anti][:,0]/(aggregated[anti][:,2]) for i,anti in enumerate(antibodies)]
violin = plt.violinplot(data,showmedians=True)
for pc in violin["bodies"]:
    pc.set_facecolor("red")
    pc.set_edgecolor("black")
    pc.set_alpha(0.5)
for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    vp = violin[partname]
    vp.set_edgecolor("black")
    vp.set_linewidth(3)
plt.xticks(np.arange(1,len(antibodies)+1,1),antibodies,fontsize=50,rotation=45)
plt.yticks(fontsize=40)
plt.xlabel("Protein",fontsize=40)
plt.ylabel("Normalised intensity",fontsize=40)
#plt.ylim([-0.1,1.5])

################### Cell quantification ##############################


### MSNs 
quant = pd.read_csv("cell_quantification.csv")

antibodies = np.array(quant.iloc[:,1:4])
antibodies = np.unique(antibodies); antibodies = antibodies[antibodies!="0"]
## aggregate intensities by receptor
channels = ["a488","a555","a647"]

aggregated = {}
plt.figure(figsize=(20,15))
for i,antibody in enumerate(antibodies):
    aggregated[antibody] = np.zeros((1,1))
    for j,channel in enumerate(channels):
        array = quant[quant[channel]==antibody]
        if len(array)!= 0:
            aggregated[antibody] = np.append(aggregated[antibody],array.iloc[:,j+7])
    aggregated[antibody] = np.delete(aggregated[antibody],0,0)
    plt.scatter(jitter(np.repeat(i+1,len(aggregated[antibody]))),aggregated[antibody],color="black",s=60)
    #plt.bar(i,np.mean(aggregated[antibody]),yerr=np.std(aggregated[antibody]),color="red",alpha=0.5)
data = [aggregated[anti] for i,anti in enumerate(antibodies)]
violin = plt.violinplot(data,showmedians=True)
for pc in violin["bodies"]:
    pc.set_facecolor("red")
    pc.set_edgecolor("black")
    pc.set_alpha(0.5)
for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    vp = violin[partname]
    vp.set_edgecolor("black")
    vp.set_linewidth(3)
    
plt.xticks(np.arange(1,i+2,1),antibodies,fontsize=40)
plt.yticks(fontsize=40)
plt.xlabel("Protein",fontsize=40)
plt.ylabel("Proportion of expression \n (Σpositive/Σnuclei)",fontsize=40)
###########################################

####### DANs#################
quant = pd.read_csv("quantification_DA.csv")

antibodies = np.array(quant.iloc[:,1:4])
antibodies = np.unique(antibodies); antibodies = antibodies[antibodies!="0"]
## aggregate intensities by receptor
channels = ["a488","a555","a594"]

aggregated = {}
plt.figure(figsize=(20,15))
for i,antibody in enumerate(antibodies):
    aggregated[antibody] = np.zeros((1,1))
    for j,channel in enumerate(channels):
        array = quant[quant[channel]==antibody]
        if len(array)!= 0:
            aggregated[antibody] = np.append(aggregated[antibody],array.iloc[:,j+7])
    aggregated[antibody] = np.delete(aggregated[antibody],0,0)
    plt.scatter(jitter(np.repeat(i+1,len(aggregated[antibody]))),aggregated[antibody],color="black",s=60)
    #plt.bar(i,np.mean(aggregated[antibody]),yerr=np.std(aggregated[antibody]),color="red",alpha=0.5)
data = [aggregated[anti] for i,anti in enumerate(antibodies)]
violin = plt.violinplot(data,showmedians=True)
for pc in violin["bodies"]:
    pc.set_facecolor("red")
    pc.set_edgecolor("black")
    pc.set_alpha(0.5)
for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
    vp = violin[partname]
    vp.set_edgecolor("black")
    vp.set_linewidth(3)
    
plt.xticks(np.arange(1,i+2,1),antibodies,fontsize=40)
plt.yticks(fontsize=40)
plt.xlabel("Protein",fontsize=40)
plt.ylabel("Proportion of expression \n (Σpositive/Σnuclei)",fontsize=40)


from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

import itertools

combs = list(itertools.combinations(range(len(antibodies[:-1])-1),2))
pcombs = list(itertools.combinations(antibodies[:-1],2)); pcombs = [list(i) for i in pcombs]
pcombcombs = list(itertools.combinations(pcombs,2)); pcombcombs = [list(i) for i in pcombcombs]
alphas = []
for i in range(len(pcombcombs)):
    a = coloc[(coloc["Image A"]==pcombcombs[i][0][0]) & (coloc["Image B"]==pcombcombs[i][0][1]) | (coloc["Image B"]==pcombcombs[i][0][0]) & (coloc["Image A"]==pcombcombs[i][0][1])]                   
    b = coloc[(coloc["Image A"]==pcombcombs[i][1][0]) & (coloc["Image B"]==pcombcombs[i][1][1]) | (coloc["Image B"]==pcombcombs[i][1][0]) & (coloc["Image A"]==pcombcombs[i][1][1])]                   
    
    alphas.append(ttest_ind(a["Area Overlap"],b["Area Overlap"])[1])
stat_array = multipletests(alphas,method='bonferroni')[1]
bool_array = multipletests(alphas,method='bonferroni')[0]

output = np.zeros((len(pcombcombs),6)).astype("object")
output[:,0:4] = np.reshape(pcombcombs,(len(bool_array),4))
output[:,4] = stat_array
output[:,5] = bool_array

### Save statistics for appendix
np.savetxt("MSN-overlap-bonferroni.csv",output,delimiter=",",fmt='%s')

for i, anti in enumerate(antibodies):    
    aggregated[anti] = np.zeros((1,4))
    for j in range(channels):
        rec_dat = quant.loc[quant.iloc[:,1+j] == anti]
        if len(rec_dat)!= 0:
            rec_dat = rec_dat.iloc[:,[i]
            aggregated[anti] = np.append(aggregated[anti],rec_dat,0)
            aggregated[anti] = aggregated[anti][1:,:]
    
    normalised = aggregated[anti][:,0]/(aggregated[anti][:,2]-aggregated[anti][:,3])
    x = jitter(np.repeat(i+1,len(aggregated[anti])),0.1)
    plt.errorbar(x, normalised, yerr=aggregated[anti][:,1]/(aggregated[anti][:,2]-aggregated[anti][:,3]), fmt="o",c="black")
    plt.scatter(x, normalised,c="black",s=40)


### Neurite plotting 
imdir = "/home/caiusgibeily/Documents/Oxdrive/Project_2/Im_analysis/MSNs/Neurite"
os.chdir(imdir)
neurites = pd.read_csv("MSN_neurites.csv",sep=",")


plt.figure(figsize=(20,15))
markers = np.append(antibodies,"DAPI")
markers = markers[markers!="DARPP32"]
aggregated = {}
for i,marker in enumerate(markers):
    aggregated[marker] = np.zeros((1,1))

    array = neurites[neurites["Antibody"]==marker]
    aggregated[marker] = np.append(aggregated[marker],array.iloc[:,3])
    aggregated[marker] = np.delete(aggregated[marker],0,0)
    plt.scatter(jitter(np.repeat(i,len(aggregated[marker]))),aggregated[marker],color="black",s=80)
    plt.bar(i,np.median(aggregated[marker]),yerr=np.std(aggregated[marker]),color="red",alpha=0.5,error_kw=dict(lw=5, capsize=5, capthick=3))
plt.xticks(np.arange(0,i+1,1),markers,fontsize=40)
plt.yticks(fontsize=40)
plt.xlabel("Marker",fontsize=40)
plt.ylabel("Proportion neurite signal ($F_{neurite}$/$F_{total}$)",fontsize=40)
###
###################### Colocalisation analysis###############################
######################## MSNs ####################################

coloc = pd.read_csv("colocalisation.csv")

##MSN_1##
imdir = "/home/caiusgibeily/Documents/Oxdrive/Project_2/Im_analysis/"
os.chdir(imdir)
FOVs = 25

red = np.array([["D1","D1","D1","M4"],["D1","D1","D1","M4"],["D1","D1","D1","M4"],["D2","D2","M4","0"],["D2","D2","M4","0"],["D2","D2","M4","0"]])
green = np.array([["a","M4","CHRNB2","OPRM1"],["a","M4","CHRNB2","OPRM1"],["a","M4","CHRNB2","OPRM1"],["D1","OPRM1","CHRNB2","0"],["D1","OPRM1","CHRNB2","0"],["D1","OPRM1","CHRNB2","0"]])
blue = np.array([["DAPI","DAPI","DAPI","DAPI"],["DAPI","DAPI","DAPI","DAPI"],["DAPI","DAPI","DAPI","DAPI"],["DAPI","DAPI","DAPI","0"],["DAPI","DAPI","DAPI","0"],["DAPI","DAPI","DAPI","0"]])

### Specify antibody correspondence - iteration from harmony is numeric: E2, F2, G2, E3, F3, G3 etc.

red_bywell = np.ndarray.flatten(red.T); red_bywell = red_bywell[red_bywell!="0"]
red_ind = np.repeat(red_bywell,FOVs)
#
green_bywell = np.ndarray.flatten(green.T); green_bywell = green_bywell[green_bywell!="0"]
green_ind = np.repeat(green_bywell,FOVs)
###
blue_ind = np.repeat("DAPI",525)


# MSN data
## MSN1
coloc.iloc[0:len(red_ind),1] = red_ind
coloc.iloc[0:len(red_ind),2] = green_ind
coloc.iloc[525:525+len(red_ind),1] = red_ind
coloc.iloc[525:525+len(red_ind),2] = blue_ind
coloc.iloc[525*2:525*2+len(red_ind),1] = green_ind
coloc.iloc[525*2:525*2+len(red_ind),2] = blue_ind
## MSN2
red = np.array([["D1","D1","D2","M4","CHRNB2"],["D1","D1","D2","M4","CHRNB2"],["D1","D1","a","M4","CHRNB2"],["D1","D1","a","M4","CHRNB2"],["D1","D2","a","M4","OPRM1"],["D1","D2","a","M4","OPRM1"]])
green = np.array([["a","OPRM1","CHRNB2","CHRNB2","OPRM1"],["a","OPRM1","CHRNB2","CHRNB2","OPRM1"],["M4","DARPP32","OPRM1","OPRM1","DARPP32"],["M4","DARPP32","OPRM1","OPRM1","DARPP32"],["CHRNB2","D1","DARPP32","DARPP32","DARPP32"],["CHRNB2","D1","DARPP32","DARPP32","DARPP32"]])


red_bywell = np.ndarray.flatten(red.T); red_bywell = red_bywell[red_bywell!="0"]
red_ind = np.repeat(red_bywell,FOVs)
#
green_bywell = np.ndarray.flatten(green.T); green_bywell = green_bywell[green_bywell!="0"]
green_ind = np.repeat(green_bywell,FOVs)
###
blue_ind = np.repeat("DAPI",750)
yellow_ind = np.repeat("CTIP2",750)

antibodies = np.unique(np.append(red,green)); antibodies = antibodies[antibodies!="0"]; antibodies = np.insert(antibodies,-1,"CTIP2");antibodies = np.insert(antibodies,-1,"DAPI")


lena = 525
lent = len(red_ind)
coloc.iloc[lena*3:lena*3+lent,1] = red_ind
coloc.iloc[lena*3:lena*3+lent,2] = green_ind
coloc.iloc[lena*3+lent:lena*3+lent+lent,1] = red_ind
coloc.iloc[lena*3+lent:lena*3+lent+lent,2] = blue_ind
coloc.iloc[lena*3+lent*2:lena*3+lent*2+lent,1] = green_ind
coloc.iloc[lena*3+lent*2:lena*3+lent*2+lent,2] = blue_ind
coloc.iloc[lena*3+lent*3:lena*3+lent*3+lent-1,1] = red_ind[:-1]
coloc.iloc[lena*3+lent*3:lena*3+lent*3+lent-1,2] = yellow_ind[:-1]
coloc.iloc[lena*3+lent*3+lent-1:lena*3+lent*3+lent-1+lent-1,1] = green_ind[:-1]
coloc.iloc[lena*3+lent*3+lent-1:lena*3+lent*3+lent-1+lent-1,2] = yellow_ind[:-1]
coloc.iloc[lena*3+lent*4+lent-2:lena*3+lent*4+lent-2+lent-1,1] = yellow_ind[:-1]
coloc.iloc[lena*3+lent*4+lent-2:lena*3+lent*4+lent-2+lent-1,2] = blue_ind[:-1]

coloc_MSN = coloc.copy()


### Heatplot
import seaborn as sns   
param = "Pearson's Coefficient"
result = coloc.groupby(["Image A", "Image B"])[param].mean().reset_index().pivot_table(index="Image A", columns="Image B", values=param,fill_value=-0.1)

result = result.combine_first(result.T)
#resultstd = coloc.groupby(["Image A", "Image B"])["Pearson's Coefficient"].ste().reset_index().pivot_table(index="Image A", columns="Image B", values="Pearson's Coefficient",fill_value=-0.1)
## Symmetric layoutttest_ind(a["Pearson's Coefficient"],b["Pearson's Coefficient"])


for i,protein in enumerate(list(result)):
    for j,protein in enumerate(list(result)):
        if result.iloc[i,j] == -0.1:
            result.iloc[i,j] = result.iloc[j,i]
np.fill_diagonal(result.to_numpy(), 0)
       
title = param
plt.figure(figsize=(18,15))
cmap = sns.color_palette("mako", as_cmap=True)
heat = sns.heatmap(result.iloc[:-1,:-1],vmin=0,vmax=1,cmap=cmap)
plt.xticks(fontsize=40,rotation=45)
plt.yticks(fontsize=40,rotation=0)
cbar = heat.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=30)
heat.collections[0].colorbar.set_label(title,fontsize=30)
plt.xlabel("Marker B",fontsize=40)
plt.ylabel("Marker A",fontsize=40)
plt.scatter(coloc["Thresholded M1"],coloc["Area Overlap"])
result_MSN = result.copy()
#####################################
##################### DANs###########

imdir = "/home/caiusgibeily/Documents/Oxdrive/Project_2/Im_analysis/"
os.chdir(imdir)
coloc_DA = pd.read_csv("colocalisation_DANs.csv")

FOVs = 25

green = np.array([["D1","D1","D1","D1","a","a"],["D2","D2","D2","D2","a","a"],["M4","M4","M4","M4","0","0"],["CHRNB2","CHRNB2","D1","D1","0","0"],["0","0","a","0","0","0"],["D1","D1","0","0","0","0"]])
red = np.array([["D2","D2","CHRNB2","CHRNB2","a","a"],["M4","M4","CHRNB2","CHRNB2","a","a"],["CHRNB2","CHRNB2","OPRM1","OPRM1","0","0"],["OPRM1","OPRM1","M4","M4","0","0"],["0","0","a","0","0","0"],["OPRM1","OPRM1","0","0","0","0"]])
blue = np.array([["DAPI","DAPI","DAPI","DAPI","a","a"],["DAPI","DAPI","DAPI","DAPI","a","a"],["DAPI","DAPI","DAPI","DAPI","0","0"],["DAPI","DAPI","DAPI","DAPI","0","0"],["0","0","a","0","0","0"],["DAPI","DAPI","0","0","0","0"]])
yellow = np.array([["TH","TH","TH","TH","a","a"],["TH","TH","TH","TH","a","a"],["TH","TH","TH","TH","0","0"],["TH","TH","TH","TH","0","0"],["0","0","a","0","0","0"],["TH","TH","0","0","0","0"]])

### Specify antibody correspondence - iteration from harmony is numeric: E2, F2, G2, E3, F3, G3 etc.

red_bywell = np.ndarray.flatten(red.T); red_bywell = red_bywell[red_bywell!="0"]
red_ind = np.repeat(red_bywell,FOVs)
#
green_bywell = np.ndarray.flatten(green.T); green_bywell = green_bywell[green_bywell!="0"]
green_ind = np.repeat(green_bywell,FOVs)
###
blue_bywell = np.ndarray.flatten(blue.T); blue_bywell = blue_bywell[blue_bywell!="0"]
blue_ind = np.repeat(blue_bywell,FOVs)
###
yellow_bywell = np.ndarray.flatten(yellow.T); yellow_bywell = yellow_bywell[yellow_bywell!="0"]
yellow_ind = np.repeat(yellow_bywell,FOVs)

antibodies = np.unique(np.append(red,green)); antibodies = antibodies[antibodies!="0"]; antibodies = np.insert(antibodies,-1,"TH");antibodies = np.insert(antibodies,-1,"DAPI")

## DA
lent = len(red_ind)
coloc_DA.iloc[0:lent,1] = red_ind
coloc_DA.iloc[0:lent,2] = green_ind
coloc_DA.iloc[lent:lent+lent,1] = red_ind
coloc_DA.iloc[lent:lent+lent,2] = blue_ind
coloc_DA.iloc[lent*2:lent*2+lent,1] = green_ind
coloc_DA.iloc[lent*2:lent*2+lent,2] = blue_ind
coloc_DA.iloc[lent*3:lent*3+lent,1] = red_ind
coloc_DA.iloc[lent*3:lent*3+lent,2] = yellow_ind
coloc_DA.iloc[lent*4:lent*4+lent,1] = green_ind
coloc_DA.iloc[lent*4:lent*4+lent,2] = yellow_ind
coloc_DA.iloc[lent*5:lent*5+lent,1] = yellow_ind
coloc_DA.iloc[lent*5:lent*5+lent,2] = blue_ind

### Heatplot

param = "Pearson's Coefficient"
result = coloc_DA.groupby(["Image A", "Image B"])[param].mean().reset_index().pivot_table(index="Image A", columns="Image B", values=param,fill_value=-0.1)

result = result.combine_first(result.T)
#resultstd = coloc.groupby(["Image A", "Image B"])["Pearson's Coefficient"].ste().reset_index().pivot_table(index="Image A", columns="Image B", values="Pearson's Coefficient",fill_value=-0.1)
## Symmetric layout


for i,protein in enumerate(list(result)):
    for j,protein in enumerate(list(result)):
        if i != 3:
            if result.iloc[i,j] == -0.1:
                result.iloc[i,j] = result.iloc[j,i]
np.fill_diagonal(result.to_numpy(), 0)
result_DA = result.copy()
title = "Overlap area (px²)"
plt.figure(figsize=(18,15))
cmap = sns.color_palette("mako", as_cmap=True)
heat = sns.heatmap(result.iloc[:-1,:-1],vmin=0,vmax=1,cmap=cmap)
plt.xticks(fontsize=40,rotation=45)
plt.yticks(fontsize=40)
cbar = heat.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=30)
heat.collections[0].colorbar.set_label(title,fontsize=30)
plt.title("Overlap area",fontsize=40)
plt.xlabel("Marker B",fontsize=40)
plt.ylabel("Marker A",fontsize=40)


from sklearn.datasets import make_regression
from sklearn.linear_model import HuberRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from matplotlib import pyplot
 
combs = list(itertools.combinations(range(len(antibodies[:-1])-1),2))
pcombs = list(itertools.combinations(antibodies[:-1],2)); pcombs = [list(i) for i in pcombs]
pcombcombs = list(itertools.combinations(pcombs,2)); pcombcombs = [list(i) for i in pcombcombs]
alphas = []

coloc = coloc_MSN ## Set for DA 
for i in range(len(pcombcombs)):
    a = coloc[((coloc["Image A"]==pcombcombs[i][0][0]) & (coloc["Image B"]==pcombcombs[i][0][1])) | ((coloc["Image B"]==pcombcombs[i][0][0]) & (coloc["Image A"]==pcombcombs[i][0][1]))]                   
    b = coloc[((coloc["Image A"]==pcombcombs[i][1][0]) & (coloc["Image B"]==pcombcombs[i][1][1])) | ((coloc["Image B"]==pcombcombs[i][1][0]) & (coloc["Image A"]==pcombcombs[i][1][1]))]                   
    
    alphas.append(ttest_ind(a["Pearson's Coefficient"],b["Pearson's Coefficient"])[1])
stat_array = multipletests(alphas,method='bonferroni')[1]
bool_array = multipletests(alphas,method='bonferroni')[0]
thresh = multipletests(alphas,method='bonferroni')[3]
output = np.zeros((len(pcombcombs),6)).astype("object")
output[:,0:4] = np.reshape(pcombcombs,(len(bool_array),4))
output[:,4] = stat_array
output[:,5] = bool_array

### Visualise output
out = pd.DataFrame(output)
out["Comb A"] = out[0]+" "+out[1]
out["Comb B"] = out[2]+" "+out[3]
out[4][np.isnan(np.array(out[4]).astype("float"))] = 1
out[6] = np.log10(np.array(out[4]).astype("float"))

result = out.groupby(["Comb A","Comb B"])[6].mean().reset_index().pivot_table(index="Comb A", columns="Comb B", values=6,fill_value=-0.1)
result_MSN = result.combine_first(result.T)


for i,stat in enumerate(list(result)):
    for j,value in enumerate(list(result)):
        if result.iloc[i,j] == -0.1:
            result.iloc[i,j] = result.iloc[j,i]
np.fill_diagonal(result.to_numpy(), 1)

title = "P-values (corrected)"
plt.figure(figsize=(18,15))
cmap=LinearSegmentedColormap.from_list('rg',["r", "w", "g"], N=256) 
norm = TwoSlopeNorm(vcenter=np.log10(thresh),vmax=1)
heat = sns.clustermap(result.iloc[:-1,:-1],vmax=1,cmap=cmap,norm=norm,dendrogram_ratio=(.2, .2),cbar_pos=(0, .2, .03, .4),annot_kws={"size": 30})

plt.xticks(fontsize=40,rotation=90)
#plt.yticks(fontsize=40)
#cbar = heat.collections[0].colorbar
heat.ax_cbar.set_ylabel("log10-P-value",size=15)
heat.ax_cbar.set_title("")
#heat.ax_cbar.set_yticks(ticks=np.arange(1,-30,-3))
# here set the labelsize by 20
#heat.collections[0].colorbar.set_label(title,fontsize=30)
#plt.title("Overlap area",fontsize=40)
#plt.xlabel("Combination B",fontsize=40)
#plt.ylabel("Combination A",fontsize=40)
heat.ax_heatmap.set_xlabel("Combination A",fontsize=20)
heat.ax_heatmap.set_ylabel("Combination B",fontsize=20)

#heat.ax_cbar.set_position((0.8, .4, .03, .4))
heat.ax_row_dendrogram.set_visible(False)
heat.ax_col_dendrogram.set_visible(False)
## Save as an svg
plt.savefig("MSN_Pearsons.svg",format="svg")

### Save statistics for appendix
np.savetxt("DA-M2-bonferroni.csv",output,delimiter=",",fmt='%s')
## Comparison 
shared = ["D1","D2","OPRM1","CHRNB2","M4","DAPI"]


share_DA = result_DA[shared]
share_DA = share_DA.loc[shared]

share_MSN = result_MSN[shared]
share_MSN = share_MSN.loc[shared]

share_MSN = share_MSN-share_MSN["DAPI"]
share_DA = share_DA-share_DA["DAPI"]


diff = np.subtract(share_MSN,share_DA)

plt.figure(figsize=(18,15))
cmap=LinearSegmentedColormap.from_list('rg',["r", "w", "g"], N=256) 
norm = TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
np.fill_diagonal(diff.to_numpy(), -100)

heat = sns.heatmap(diff.iloc[:-1,:-1],cmap=cmap,norm=norm)
cbar = heat.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=30)
plt.xlabel("Combination B",fontsize=40)
plt.ylabel("Combination A",fontsize=40)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

plt.title("Pairwise correlation difference (MSN-DAN)",fontsize=30)
## Regression of area/correlation coefficients


correlation_reg = plt.scatter(coloc_DA["Thresholded M1"],coloc_DA["Pearson's Coefficient"])

