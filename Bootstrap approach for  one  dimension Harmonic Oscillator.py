#!/usr/bin/env python
# coding: utf-8

# In[1]:


# bootstrap for one dimension harmonic oscillator.
# here we are trying to plot tanh(detM) vs energy and see that at higher  order of matrix the accuracy is good.
import numpy as np
import matplotlib.pyplot as plt

def rec(s, E):
    if s == 0:
        return 1
    elif s == 2:
        return E
    elif s % 2 == 1:
        return 0
    else:
        return 2 * E * (s - 1) * rec(s - 2, E) / s + ((s - 1) * (s - 2) * (s - 3)*rec(s-4,E)) / (4 * s)

def create_matrix(n, E):
    """
    Create the matrix M with elements M_ij = <x^(i+j)>.

    Args:
    - n (int): Size of the matrix (n x n).
    - E (float): Value of E.

    Returns:
    - numpy.ndarray: The matrix M.
    """
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = rec(i+j,E)
    return M

n=int(input("Enter the order of matrix:"))
d=int(input("Enter range of E:"))

# Define the step size
step_size = 0.2

# Initialize arrays to store E and tanh_det_M values
E_values = []
tanh_det_M_values = []

# Loop over the values of E with a step size of 0.2
for E in np.arange(0, d , step_size):
   
    M = create_matrix(n, E)
    det_M = np.linalg.det(M)
    tanh_det_M = np.tanh(det_M)
    
    E_values.append(E)
    tanh_det_M_values.append(tanh_det_M)
    print(tanh_det_M, E)

# Plot E against tanh_det_M
plt.plot(E_values, tanh_det_M_values)
plt.xlabel('E')
plt.ylabel('tanh(det(M))')
plt.title('Plot of tanh(det(M)) vs E')
plt.grid(True)
plt.show()


# In[ ]:




