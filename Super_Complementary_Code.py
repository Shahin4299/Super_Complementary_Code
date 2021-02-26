#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import copy as cp
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import correlate
# from scipy import signal
import os


# ### Import Walsh Hadamard code as our initial value of matrices

# In[2]:


def WalshGenerator(code_length):
        
        code_length
        code=[[-1,-1],[-1,1]]
        [r1,c1]=np.shape(code)
        while r1<code_length:
            code=np.concatenate((code,code))
            code=np.concatenate((code,code),axis=1)
            for i in range(r1):
                for j in range(c1):
                    code[i+r1][j+c1]=-1*code[i][j]
            [r1, c1] = np.shape(code)
        #print(code)
        return code


# ### Specify the dimension of initial values and calculate C from a & b matrices - multiply the whole row of a to matrix b element

# In[3]:


dim = 2
a=np.array(WalshGenerator(dim))
b=np.array(WalshGenerator(dim))
d=np.array(WalshGenerator(dim))

n1,n2=np.shape(a)
c=np.array([])
e=np.array([])
t=np.array([])
s=np.array([])
eext=np.array([])
for i in range(n1):
    for j in range(n1):
        c=np.append(c,b[j][i]*a[j][:])
c=np.reshape(c,(n1,n1**2))


# In[4]:


# c.shape


# ### To generate E which is Complete Complementary Code, that multiply C and our initial value D- multiply each element of C to the elements of D

# In[5]:


for i in range(n1):
    for j in range(n1):
        dd=d[j][:]
        dd = np.tile(dd,n1)
        cc = c[i][:]
        mp=np.multiply(cc,dd)
        e=np.append(e,mp)
e.shape


# ### Generating n**2 matrix from the previous equations- Complete Complementary Code generated

# In[6]:


e=np.reshape(e,(n1**2,n1**2))
# e.shape


# ### Creating Super Complementary matrix shape

# In[7]:


t=cp.copy(e)
m1,m2=np.shape(t)
td=np.array([])
s=np.empty((2*n1**2,2*n1**2))
s.shape


# ### Generate Super Complementary codes by using Complete Complementary Codes
# 

# In[8]:


for i in range(0,n1**2,2):
    td=np.vstack([t[i][:],t[i+1][:]])
    td=np.reshape(td,(1,2*n1**2), order='F')
    s[i*2]=td
    
    for j in range(2*n1**2):
        s[2*i+1][j]=(-1)**(j)*s[2*i][j]
    for k in range(0,2*n1**2-1,2):
        c1=s[2*i][k]
        c2=s[2*i][k+1]
        s[2*i+2][k]=c2
        s[2*i+2][k + 1]=c1
    for l in range(2*n1**2):
         s[2*i+3][l]=(-1)**(l)*s[2*i+2][l]

s
            


# ### Checking orthogonality of the rows

# In[9]:


if np.array([[(s[i]*s[j]).sum() if i!=j else 0 for i in range(s.shape[0])] for j in range(s.shape[0])]).sum() == 0:
    print("All rows are orthogonal")
else:
    print("Rows are Not orthogonal")


# ### Define a method for autocorrelation

# In[10]:


def autoCorr(x, lag):
    result = np.correlate(x, x, mode= 'full')
    return result[result.size//2:]/np.correlate(x,x)
#     return pd.Series(sm.tsa.acf(x, lag))


# ### Define a method for cross correlation

# In[11]:


def crossCorr(row1, row2):
    result =  np.correlate(row1, row2 , 'full')
    rx = np.correlate(row1, row1)
    ry = np.correlate(row2, row2)
    result = result/((rx*ry)**0.5)
    return result


# In[ ]:





# ### Plot the results of Auto-Corr and Cross-Corr for Super Complementary code
# 

# In[12]:


plt.figure(figsize = (50,50))
nrow_plot, ncol_plot =s.shape
lag = ncol_plot
matrix = s
for j in range(nrow_plot):
    for k in range(ncol_plot):
        plt.subplot(nrow_plot, ncol_plot, 1+j*nrow_plot+k)
        if j == k:
            ac = autoCorr(matrix[j], lag)
#             ac = mmscaler.fit_transform(np.expand_dims(ac, axis=1))
            plt.plot(ac, 'r')
            plt.plot([0,nrow_plot-1], [0,0])
        else:
            cc = crossCorr(matrix[j], matrix[k])
#             cc = mmscaler.fit_transform(np.expand_dims(cc, axis=1))
            plt.plot(cc)
            plt.plot([0,2*nrow_plot-2], [0,0])
        plt.ylim(-1,1)
plt.savefig('Correlation_Super_8_nu')


# In[13]:


cc


# ### Define a method for adding noise 

# In[14]:


clean_signal = pd.DataFrame(s, dtype=float) 
# print(clean_signal)
mu, sigma = 0, 0.1 
noise = np.random.normal(mu, sigma, s.shape)
# print(noise)
signal = clean_signal + noise
# print(signal)


# In[15]:


signal


# ### Plot the results by adding noise 

# In[16]:


plt.figure(figsize = (50,50))
nrow_plot, ncol_plot =s.shape
lag = ncol_plot
matrix = signal.values
for j in range(nrow_plot):
    for k in range(ncol_plot):
        plt.subplot(nrow_plot, ncol_plot, 1+j*nrow_plot+k)
        if j == k:
            ac = autoCorr(matrix[j], lag)
#             ac = mmscaler.fit_transform(np.expand_dims(ac, axis=1))
            plt.plot(ac, 'r')
            plt.plot([0,nrow_plot-1], [0,0])
        else:
            cc = crossCorr(matrix[j], matrix[k])
#             cc = mmscaler.fit_transform(np.expand_dims(cc, axis=1))
            plt.plot(cc)
            plt.plot([0,2*nrow_plot-2], [0,0])
        plt.ylim(-1,1)
plt.savefig('Correlation_Super_32_nu_Noise')


# ### Define a method for flipup

# In[17]:


def flipUp(mtr):
    flip = []
    rows = len(mtr)
    j = -1
    for i in range (rows):
        flip.append(mtr[j])
        #print(mtr[j])
        j -= 1
    return(flip)


# In[18]:


fliped = np.array(flipUp(s))
# fliped


# ### Define a method for randomize

# In[19]:


import random
def randomize(mtrx, num_rows, order_of_codes):
    
    if order_of_codes < 3:
        print('The number of columns is lower than expected !!!')
        return 0
    row , col = np.shape(mtrx)
    sel_rows = random.sample(range(0, row), num_rows)
    sel_cols = 2**(order_of_codes-1) # number of chips (code length)
    
    sel_mrtx = []
    for i in sel_rows:
        sel_mrtx.append(mtrx[i])
    sel_mrtx = np.array(sel_mrtx)   
    final = sel_mrtx[:, :order_of_codes]
    return final


# In[20]:


randomize(s.astype(int), 3,4)


# ##### First step to create Extended Complementary Codes, by using Complete Complementary code with code order N and length N**n

# In[21]:


ne1,ne2=np.shape(e)
G=list()
for i in range(0,ne1-1,2):
    e1=cp.copy(e[i][:])
    e2 = cp.copy(e[i+1][:])
    for j in range(ne1):
        G.append(e1[j])
        G.append(e2[j])
G=np.reshape(G,(n1,n1*ne1))
G.shape


# ### Making first list to use for further calculation of Extended Complementary code

# In[22]:


E=list()
for i in range(n1):
    for j in range(n1):
        dd=d[j][:]
        #print(d[j][:])
        dd = np.tile(dd,n1**2)
        GG= G[i][:]
        mp=np.multiply(GG,dd)
        E=np.append(E,mp)
E.shape


# ### Generated based on a Complete Complementary code set of order N and length N**n−1

# In[23]:


E=np.reshape(E,(n1**2,n1**3))
#print(E)
Gnew=list()
ns1,ns2=np.shape(E)
#print(ns1,ns2)
for i in range(0,ns1-1,2):
    E1 = cp.copy(E[i][:])
    E2 = cp.copy(E[i+1][:])
    for j in range(ns2):
        Gnew.append(E1[j])
        Gnew.append(E2[j])
Gnew=np.reshape(Gnew,(n1,n1**4))
# Gnew


# ### Generating Extended Complementary code

# In[24]:


ESC=list()
for i in range(n1):
    for j in range(n1):
        dd=d[j][:]
        dd = np.tile(dd,n1**3)
        Gnews= Gnew[i][:]
        mp=np.multiply(Gnews,dd)
        ESC=np.append(ESC,mp)
ESC=np.reshape(ESC,(n1**2,n1**4))* (-1)
ESC.shape


# ##### By controlling the extended complementary code set, we can generate a super complementary code with different element code lengths.

# In[25]:


s_ext=np.zeros([2*n1**2,2*n1**4])
for i in range(0,n1**2,2):
    for j in range(0,n1**4,n1**2):
        s_ext[i*2][2*j:2*j+n1**2]=ESC[i][j:j+n1**2]
        s_ext[i*2][2*j+n1**2:2*(j+n1**2)] =(ESC[i+1][j: j+n1**2])

        #print(2*j+n1**2)
        #print(2*(j+n1**2))
        #print(int(j/4))
        #print((int(j/4)) + n1**2)

        s_ext[i * 2+1][2 * j:2 * j + n1 ** 2] = ESC[i][j:j+n1**2]
        s_ext[i * 2+1][2 * j + n1 ** 2:2 * (j + n1 ** 2)] = (ESC[i+1][j:j+n1**2])*(-1)
        s_ext[i * 2+2][2 * j:2 * j + n1 ** 2] = ESC[i+1][j:j+n1**2]
        s_ext[i * 2+2][2 * j + n1 ** 2:2 * (j + n1 ** 2)] =ESC[i][j:j+n1**2]
        s_ext[i * 2+3][2 * j:2 * j + n1 ** 2] = ESC[i + 1][j:j+ n1 ** 2]
        s_ext[i * 2+3][2 * j + n1 ** 2:2 * (j + n1 ** 2)] = ESC[i][j:j+n1 ** 2]*(-1)

s_ext.shape


# ### Checking the orthogonality of th codes

# In[34]:


def isOrthogonal(g_array):
    total_sum_G = 0
    for i in range(len(g_array-1)):
        for j in range(i+1 , len(g_array-1)):
            total_sum_G += sumarr(g_array[i] , g_array[j])
    return total_sum_G == 0
if isOrthogonal:
    print("All rows are orthogonal")
else:
    print("Rows are not orthogonal")


# In[35]:


plt.figure(figsize = (50,50))
nrow_plot, ncol_plot =s_ext.shape
lag = ncol_plot
matrix = s_ext
for j in range(nrow_plot):
    for k in range(nrow_plot):
        plt.subplot(nrow_plot, nrow_plot, 1+j*nrow_plot+k)
        if j == k:
            ac = autoCorr(matrix[j], lag)
#             ac = mmscaler.fit_transform(np.expand_dims(ac, axis=1))
            plt.plot(ac, 'r')
            plt.plot([0,ncol_plot-1], [0,0])
        else:
            cc = crossCorr(matrix[j], matrix[k])
#             cc = mmscaler.fit_transform(np.expand_dims(cc, axis=1))
            plt.plot(cc)
            plt.plot([0,2*ncol_plot-2], [0,0])
        plt.ylim(-1,1)
plt.savefig('Correlation_Super_8_EXT_nu')


# ### Adding noise to the Super extended Complementary code

# In[36]:


clean_signal = pd.DataFrame(s_ext, dtype=float) 
# print(clean_signal)
mu, sigma = 0, 0.1 
noise = np.random.normal(mu, sigma, s_ext.shape)
# print(noise)
signal_ext = clean_signal + noise
signal_ext.shape


# ### Plot the Super Complementary code by adding noise

# In[37]:


plt.figure(figsize = (50,50))
nrow_plot, ncol_plot =s_ext.shape
lag = ncol_plot
matrix = signal_ext.values
for j in range(nrow_plot):
    for k in range(nrow_plot):
        plt.subplot(nrow_plot, nrow_plot, 1+j*nrow_plot+k)
        if j == k:
            ac = autoCorr(matrix[j], lag)
#             ac = mmscaler.fit_transform(np.expand_dims(ac, axis=1))
            plt.plot(ac, 'r')
            plt.plot([0,ncol_plot-1], [0,0])
        else:
            cc = crossCorr(matrix[j], matrix[k])
#             cc = mmscaler.fit_transform(np.expand_dims(cc, axis=1))
            plt.plot(cc)
            plt.plot([0,2*ncol_plot-2], [0,0])
        plt.ylim(-1,1)
plt.savefig('Correlation_Super_8_EXT_nu_Noise')

