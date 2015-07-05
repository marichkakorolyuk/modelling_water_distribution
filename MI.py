# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 10:56:42 2015

@author: mariakoroliuk
"""
# connecting libtaries
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

# cleaning screen
plt.clf()

n_iter=20000  # amount of times we are sampling the capilars
r_mean=5.42/2  # mean value for fiber radious
r_div=2.7/2 # this is a parameters that represent how far fiber radius deviates  from mean value
min_value_r_f=0.1  # fiber can not go smaller that this  (nanaometers)
po_dev=10  #deviation for the density of material kg/m^3
min_value_pho_b=20  # minimal value of density bottom layer
min_value_pho_t=30   # minimal value of density top layer

T=0.0728   # this is the constant at fixed temperature
g=9.8  #gravity constant m/s^2
gamma=0.1047  # contact angle translates to radians (6 degrees)
pho_water=1000 # density of water (or other liquid) kg/m^3
pho_f=2700    # density of fiber
h_max=0.12  #  this is higth of whole fiber
ratio=1  #(could be 0.75 to make it 2 layers).1 to look at only one layer
h_b=ratio*h_max   # this is  bottom layer  higth
h_t=(1-ratio)*h_max    # this is top layer   higth
pho_b_mean=40
# next 4 lines are for plotting
from matplotlib.patches import Rectangle
figure(1)


FSC=0.3
r_sc=150*10**(-6)

# main function
# if takes as argument values of top layer and bottom layer densities. If pho_t_mean is not specified , if will avtomatically be the same as pho_b_mean
# so you can run it as  SaturateWater(50) or SaturateWater(40,50) or whatever. Returns theta and h_0 (saturation and heigth)

def SaturateWater(pho_b_mean,pho_t_mean=pho_b_mean):
    h_vector=[]
    pho_b_mean=pho_b_mean
    pho_t_mean=pho_t_mean

    for i in range(n_iter): 
        r_f=(np.random.normal(r_mean, r_div))   # GENERATE fiber radius
        while r_f<min_value_r_f:
            r_f=(np.random.normal(r_mean, r_div))
        pho_b_m=(np.random.normal(pho_b_mean,po_dev))  #GENERATE fiber density!!! do not allow it to go too small
  
        while pho_b_m<min_value_pho_b:
            pho_b_m=(np.random.normal(pho_b_mean,po_dev))  #GENERATE fiber density!!!
            
        eps=1-pho_b_m/pho_f  # poriousity
        h_star=2*np.cos(gamma)*T/(pho_water*g*r_f)  # dimentionless parametes h_star
        h=h_star*(1/eps-1)

        
    # next loop is never used now with one layer!
    # next loop is for water that asceessed first layer
        if h*10**6>h_b:
            r_f=(np.random.normal(r_mean, r_div))   
            while r_f<min_value_r_f:
                r_f=(np.random.normal(r_mean, r_div))
            pho_b_m=(np.random.normal(pho_b_mean,po_dev))  
  
            while pho_b_m<min_value_pho_b:
                pho_b_m=(np.random.normal(pho_b_mean,po_dev)) 
            
            eps=1-pho_b_m/pho_f
            h_star=2*np.cos(gamma)*T/(pho_water*g*r_f)
            h=h_star*(1/eps-1)

            h=h*10**6+h_b # pressure head for a if condition (if it got to a second layer)
            
        else:h=h*10**6 #pressure head if it never got for an if condition (we are consdering one layer of water never get to second layer)

        h_vector.append(h) # 
    theta=[]  # create vector empty
    h_0=np.arange(0,h_max,0.0001) # min, max, step
    h_vector=np.array(h_vector) # change a type of vector to np.vector, that makes it similar to matlab 
    
    for h1 in h_0: 
        if h1>h_b:   # calculating sat  and residual content for current layer (top of bottom)
            theta_sat=1-pho_t_mean/pho_f   # theta saturated
            theta_res1=1.5*(float(pho_t_mean)/pho_f)*FSC*r_mean*10**(-6)*np.cos(gamma)*T/(r_sc**3*pho_water*g)

        else:
            theta_sat=1-pho_b_mean/pho_f   # theta saturated
            theta_res1=1.5*(float(pho_b_mean)/pho_f)*FSC*r_mean*10**(-6)*np.cos(gamma)*T/(r_sc**3*pho_water*g)


        theta.append((sum(h_vector>h1)/float(len(h_vector)))*(theta_sat-theta_res1)+theta_res1) # this is  probability already normilized
    
    return(theta,h_0)

theta,h_0=SaturateWater(pho_b_mean)  

#plotting
someX, someY = n_iter/2, h_b/2       
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((someX - n_iter/2, someY-h_b/2), n_iter, h_b, alpha=0.3,facecolor='y'))
someX, someY = n_iter/2, h_t/2+h_b
currentAxis.add_patch(Rectangle((someX - n_iter/2, someY-h_t/2), n_iter, h_t, alpha=0.3,facecolor='r'))


# van Genuhten

theoretical=[]
theta_sat=1-pho_b_mean/pho_f
def func(x,n,a,m,theta_res):return theta_res+(theta_sat-theta_res)/(1+(a*x**n)**(m)) # this function is returning calculation for van Genucten model. 

x0=np.array([2.7,150,1,0.2])  #start position for optimazing (could be anything alse!)

import scipy.optimize as optimization
params=optimization.curve_fit(func,h_0,theta, x0)[0] #LSF for van Genuhten.  FUNC shows what to fit, h_0 is argument, theta is function itself, x0- starting position
n=params[0] 
a=params[1]
m=params[2]
theta_res=params[3]


plot(theta,h_0,'g',linewidth=2,alpha=1,label='model')

theta=np.array(theta)  # changing vector type for a type when we could use np.mean later

for h1 in h_0:
    theoretical.append((theta_res+(theta_sat-theta_res)/(1+(a*h1**n)**(m)))) # calculating theta for van Genucten


plot(theoretical,h_0,'m',linewidth=2,alpha=1,label='van Genucten (fitted by LSM)')


s='density= '+str(pho_b_mean)+' kg/m^3'
avarage_water_content=round(np.mean(theta[h_0<h_b]),2)

s1='average WC(layer)='+str(avarage_water_content)
text(0.5,0.01, s1, fontsize=12)
text(0.05,0.02, s, fontsize=12)



avarage_water_content=round(np.mean(theta),2)

s1='average WC (total)='+str(avarage_water_content)
text(0.05,0.05, s1, fontsize=12)
xlim([0,1])


plt.xlabel('water saturation (%/100)', fontsize=14)
plt.ylabel('height (meters)', fontsize=14)    
plt.legend()
plt.show()    