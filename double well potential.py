#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
import itertools
memo={}
def calculate_x(t, E,c):
    if (t, E,c) in memo:
        return memo[(t, E,c)]
    
    if t == 0:
        return 1
    elif t == 2:
        return c
    elif t==4:
        return ((E-1.25)+2*c)*5/3
    elif t % 2 == 1:
        return 0
    else:
        term1 = 4 * (t - 3) * (E - 1.25) * calculate_x(t - 4, E,c)
        term2 = (t - 3) * (t - 4) * (t - 5) * calculate_x(t - 6, E,c)
        term3 = 4 * (t - 2) * calculate_x(t - 2, E,c)
        denominator = 4 * (t - 1)/5
        result = (term1 + term2 + term3) / denominator
        memo[(t, E,c)] = result
    return result

def create_matrix(n, E,c):
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = calculate_x(i+j, E,c)
    return M


def is_semi_pos_def(matrix):
    eigenvalues, _ = np.linalg.eig(matrix)
    return np.all(eigenvalues >= 0)
n=int(input("Enter the order of matrix:"))
d=int(input("Enter range of E:"))
c=int(input("Enter guess value of x^2:"))
range1=np.arange(0,d,0.001)
range2=np.arange(0,c,0.001)
E_Values=[]
x_values=[]


for E in range1:
    for x in range2:
        M=create_matrix(n,E,x)
        if is_semi_pos_def(M):
            E_Values.append(E)
            x_values.append(x)
E_V=[]
E_v=[]
for i in E_Values:
    if i<=1:
        E_V.append(i)
    else:
        E_v.append(i)
A=max(E_V)-min(E_V)
B=max(E_v)-min(E_v)
print("energy width at ground" ,A)
print("energy width of first excited", B)
Energy_diff=(max(E_v)+min(E_v))/2-(max(E_V)+min(E_V))/2
print("Energy difference is",":",Energy_diff)
    
plt.scatter(E_Values,x_values)
    
        


# In[ ]:




