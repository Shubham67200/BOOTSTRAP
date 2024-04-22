#!/usr/bin/env python
# coding: utf-8

# In[33]:


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

def is_positive_definite(matrix):
    try:
        cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False

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
    print(tanh_det_M," : ", E)
E_1=[]

allowed=[]    
       
for i in np.arange(0, d , step_size):
    if is_positive_definite(create_matrix(n,i)):
        E_1.append(i)
        point=1
    else:
        point=0
    allowed.append(point)
print(E_values,";",allowed)

# Plot E against tanh_det_M
plt.plot(E_values, tanh_det_M_values)
plt.xlabel('E')
plt.ylabel('tanh(det(M))')
plt.title('Plot of tanh(det(M)) vs E')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[11]:





# In[1]:


# bootstrap for one dimension harmonic oscillator.
# here we are trying to plot tanh(detM) vs energy and see that at higher  order of matrix the accuracy is good.
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
# Memoization dictionary
memo = {}


def rec(s, E):
    # Check if the result for current s and E is already calculated
    if (s, E) in memo:
        return memo[(s, E)]
    
    if s == 0:
        result = 1
    elif s == 2:
        result = E
    elif s % 2 == 1:
        result = 0
    else:
        result = 2 * E * (s - 1) * rec(s - 2, E) / s + ((s - 1) * (s - 2) * (s - 3)*rec(s-4,E)) / (4 * s)
    
    # Store the result in memoization dictionary
    memo[(s, E)] = result
    return result

def create_matrix(n, E):
    
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = rec(i+j, E)
    return M
def is_positive_definite(matrix):
    try:
        cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False



while True:
    
    n=int(input("Enter the order of matrix:"))
    d=int(input("Enter range of E:"))
        # Define the step size
    step_size = 0.00001

    # Initialize arrays to store E and tanh_det_M values
    E_values = []
    tanh_det_M_values = []
    allowed=[]

    # Loop over the values of E with a step size of 0.0005
    Range = np.arange(0, d, step_size)
    E_1=[]   
    i = 0
    while i < len(Range):
        E = Range[i]
        for j in range(1,n):
            M = create_matrix(j, E)

            E_values.append(E)

            if is_positive_definite(M):
                E_1.append(E)
                break
            else:

                Range[i] = np.nan 
        i += 1  
    E_2=[]
    for k in E_1:
        M=create_matrix(n,k)
        if is_positive_definite(M):
            point=1
            E_2.append(k)
        else:
            point=0
        allowed.append(point)


    E_3=[]
    E_4=[]
    #print(E_1)
    #print(E_2)
    for i in E_2:
        if i<=6:
            E_3.append(i)
        else:
            E_4.append(i)



    #print(Range[~np.isnan(Range)])



    print(sum(E_3))
    print(len(E_3))

    A_E=sum(E_3)/len(E_3)
    print('average energy',':',A_E)
    error=(0.5-A_E)/0.5*100
    print("percentage error is",":",error)                 

    energy_width=max(E_3)-min(E_3)  
    print("energy width"," : " ,energy_width)        





    # Plot allowed vs E
    plt.plot(E_1, allowed)
    plt.xlabel('E')
    plt.ylabel('Allowed')
    plt.title('Allowed vs E')
    plt.grid(True)

    # Add horizontal lines to indicate energy width

    # Show the energy width in the legend
      #  plt.text(min(E_1) + energy_width / 2, min(allowed) + 0.5, f'Energy Width: {energy_width:.5f}', ha='center')


    plt.show()
    choice = input("Do you want to enter another set of values? (yes/no): ")
    if choice.lower() != 'yes':
        break


# # bootstrap for anharmonic oscillator (double well potential)

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
import itertools
memo={}

def calculate_x(t, E,c):
    g=1
    V=6.25
    m=5
 
    if (t, E,c) in memo:
        return memo[(t, E,c)]
    
    if t == 0:
        return 1
    elif t == 2:
        return c
    elif t==4:
        return ((E-6.25)+2*5*c)/3
    elif t % 2 == 1:
        return 0
    else:
        term1 = 4 * (t - 3) * (E - V) * calculate_x(t - 4, E,c)
        term2 = (t - 3) * (t - 4) * (t - 5) * calculate_x(t - 6, E,c)
        term3 = 4 * 5*(t - 2) * calculate_x(t - 2, E,c)
        denominator = 4 * (t - 1)
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
while True:
    
    n=int(input("Enter the order of matrix:"))
    d=int(input("Enter range of E:"))
    c=int(input("Enter guess value of x^2:"))
    range1=np.arange(2,d,0.0001)
    range2=np.arange(1.5,c,0.0001)
    E_Values=[]
    x_values=[]


    for E in range1:
        for x in range2:
            M=create_matrix(n,E,x)
            if is_semi_pos_def(M):
                E_Values.append(E)
                x_values.append(x)
    
    plt.scatter(E_Values,x_values)
    choice = input("Do you want to enter another set of values? (yes/no): ")
    if choice.lower() != 'yes':
        break
    
        


# In[14]:


E1=[]
E2=[]
#print(E_Values)
for i in E_Values:
    if i<=3:
        E1.append(i)
    else:
        E2.append(i)
#print(E1)
#print(E2)
a=sum(E1)/len(E1)
b=sum(E2)/len(E2)
print("average energy at lower :",sum(E1)/len(E1),"+0.001")
print("average energy at upper :",sum(E2)/len(E2),"+0.001")
error = (b-a)*0.001*(1/float(a)+1/float(b))
print("splitting between energy :",b-a,"+",error)


# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:


print(25/4)


# In[5]:


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
    
plt.scatter(E_Values,x_values,marker='o')
    
        


# In[4]:


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
    
        


# In[23]:


print(sum(E_v))
print(len(E_v))
a='{:.3f}'.format(sum(E_v)/len(E_v))
print("upper energy ",":",a)


# In[24]:


print('{:.3f}'.format(sum(E_V)))
print(len(E_V))
b='{:.3f}'.format(sum(E_V)/len(E_V))
print("lower energy ",":",b)


# In[22]:


error = Energy_diff*0.001*(1/float(a)+1/float(b))
print(error)


# In[16]:


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
    
        


# In[25]:


print('{:.3f}'.format(sum(E_V)))
print(len(E_V))
b='{:.3f}'.format(sum(E_V)/len(E_V))
print("lower energy ",":",b)


# In[26]:


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
    
plt.scatter(E_Values,x_values)
plt.xlabel('E')
plt.ylabel('x^2')
plt.title('Scatter plot of E vs x^2')
plt.legend()

# Set limits for x and y axes to include full range of values
plt.xlim(min(range1), max(range1))
plt.ylim(min(range2), max(range2))

plt.show()
        


# In[27]:


print(E_Values)


# In[2]:


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
    
plt.scatter(E_Values,x_values)
plt.xlabel('E')
plt.ylabel('x^2')
plt.title('Scatter plot of E vs x^2')
plt.legend()

# Set limits for x and y axes to include full range of values
plt.xlim(min(range1), max(range1))
plt.ylim(min(range2), max(range2))

plt.show()
        


# In[4]:





# In[17]:


E_values[]


# In[37]:


len(E_new)


# In[2]:


# bootstrap for one dimension harmonic oscillator.
# here we are trying to plot tanh(detM) vs energy and see that at higher  order of matrix the accuracy is good.
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
# Memoization dictionary
memo = {}


def rec(s, E,a):
    if a==1:
        
    # Check if the result for current s and E is already calculated
        if (s, E,a) in memo:
            return memo[(s, E,a)]
    
        if s == 0:
            result = 1
        elif s == 2:
            result = E
        elif s % 2 == 1:
            result = 0
        else:
            result = 2 * E * (s - 1) * rec(s - 2, E) / s + ((s - 1) * (s - 2) * (s - 3)*rec(s-4,E)) / (4 * s)
    
    # Store the result in memoization dictionary
        memo[(s, E,a)] = result
        return result
    if a==2:
        if (s, E,a) in memo:
            return memo[(s, E,a)]
        if s==0:
            result =1
        elif s==2:
            return c
        elif s==4:
            return E/5
        else:
            term1=2*E*(s-3)*rec(s-4,E)
            term2=(s-3)*(s-4)*(s-5)*rec(s-6,E)/4
            denominator=2*(s+1)
            result =(term1+term2)/denominator
        memo[(s,E,a)]=result
        return result

        
        

def create_matrix(n, E):
    
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = rec(i+j, E)
    return M
def is_positive_definite(matrix):
    try:
        cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False




n=int(input("Enter the order of matrix:"))
d=int(input("Enter range of E:"))
    # Define the step size
step_size = 0.00001

# Initialize arrays to store E and tanh_det_M values
E_values = []
tanh_det_M_values = []
allowed=[]

# Loop over the values of E with a step size of 0.0005
Range = np.arange(0, d, step_size)
E_1=[]   
i = 0
while i < len(Range):
    E = Range[i]
    for j in range(1,n):
        M = create_matrix(j, E)
            
        E_values.append(E)
           
        if is_positive_definite(M):
            E_1.append(E)
            break
        else:
            
            Range[i] = np.nan 
    i += 1  
E_2=[]
for k in E_1:
    M=create_matrix(n,k)
    if is_positive_definite(M):
        point=1
        E_2.append(k)
    else:
        point=0
    allowed.append(point)


E_3=[]
E_4=[]
#print(E_1)
#print(E_2)
for i in E_2:
    if i<=6:
        E_3.append(i)
    else:
        E_4.append(i)
    
        
        
#print(Range[~np.isnan(Range)])


       
print(sum(E_3))
print(len(E_3))
    
A_E=sum(E_3)/len(E_3)
print('average energy',':',A_E)
error=(0.5-A_E)/0.5*100
print("percentage error is",":",error)                 
                     
energy_width=max(E_3)-min(E_3)  
print("energy width"," : " ,energy_width)        





# Plot allowed vs E
plt.plot(E_1, allowed)
plt.xlabel('E')
plt.ylabel('Allowed')
plt.title('Allowed vs E')
plt.grid(True)

# Add horizontal lines to indicate energy width

# Show the energy width in the legend
  #  plt.text(min(E_1) + energy_width / 2, min(allowed) + 0.5, f'Energy Width: {energy_width:.5f}', ha='center')


plt.show()
    


# In[1]:


def rec(m, E, a):
    if m == 0:
        return 1
    elif m == 2*a:
        return E/2
    elif m==2*(a+1):
        return 3*E*rec(2,E,a)/(a+3)+3/(4*(a+3))
    elif m % 2 == 1:
        return 0
    else:
        term1 = 2*(m-2*a+1)*E*rec(m-2*a, E, a)
        term2 = ((m-2*a)*(m-2*a-1)*(m-2*a+1)*rec(m-2*a-2, E, a))/4
        denominator = 2*(m-a+1)
        result = (term1 + term2) / denominator
    return result


# In[ ]:


rec(6,1,2)


# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
import itertools
memo = {}


def rec(s, E):
    # Check if the result for current s and E is already calculated
    if (s, E) in memo:
        return memo[(s, E)]
    
    if s == 0:
        result = 1
    elif s == 2:
        result = E/2
    elif s % 2 == 1:
        result = 0
    else:
        result =  E * (s - 1) * rec(s - 2, E) / s + ((s - 2) *( (s - 1) * (s - 3)/4-2*l*(l+1))*rec(s-4,E))/ ( 2* s)
    
    # Store the result in memoization dictionary
    memo[(s, E)] = result
    return result
def create_matrix(n, E):
    
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = rec(i+j, E)
    return M
def is_positive_definite(matrix):
    try:
        cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False




n=int(input("Enter the order of matrix:"))
d=int(input("Enter range of E:"))
l=int(input("enter the range of l:"))
    # Define the step size
step_size = 0.00001

# Initialize arrays to store E and tanh_det_M values
E_values = []
tanh_det_M_values = []
allowed=[]

# Loop over the values of E with a step size of 0.0005
Range = np.arange(0, d, step_size)
E_1=[]   
i = 0
while i < len(Range):
    E = Range[i]
    for j in range(1,n):
        M = create_matrix(j, E)
            
        E_values.append(E)
           
        if is_positive_definite(M):
            E_1.append(E)
            break
        else:
            
            Range[i] = np.nan 
    i += 1  
E_2=[]
for k in E_1:
    M=create_matrix(n,k)
    if is_positive_definite(M):
        point=1
        E_2.append(k)
    else:
        point=0
    allowed.append(point)


E_3=[]
E_4=[]
#print(E_1)
#print(E_2)
for i in E_2:
    if i<=6:
        E_3.append(i)
    else:
        E_4.append(i)
    
        
        
#print(Range[~np.isnan(Range)])


       
print(sum(E_3))
print(len(E_3))
    
#A_E=sum(E_3)/len(E_3)
#print('average energy',':',A_E)
#error=(0.5-A_E)/0.5*100
#print("percentage error is",":",error)                 
                     
energy_width=max(E_3)-min(E_3)  
print("energy width"," : " ,energy_width)        





# Plot allowed vs E
plt.plot(E_1, allowed)
plt.xlabel('E')
plt.ylabel('Allowed')
plt.title('Allowed vs E')
plt.grid(True)

# Add horizontal lines to indicate energy width

# Show the energy width in the legend
  #  plt.text(min(E_1) + energy_width / 2, min(allowed) + 0.5, f'Energy Width: {energy_width:.5f}', ha='center')


plt.show()
    


# In[ ]:


## 


# In[ ]:





# In[ ]:




