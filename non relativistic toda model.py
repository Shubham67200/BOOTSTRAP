#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import math
from multiprocessing import Pool

# Memoization dictionary
memo = {}

# Recursive function to compute values
def rec(m, n, E, h, f):
    if (m, n, E, h, f) in memo:
        return memo[(m, n, E, h, f)]
    
    if n < 0:
        result = rec(m, -n, E, h, f) / (-1) ** m
    elif m == 0:
        if n == 0:
            result = 1
        elif n == 1:
            result = f
        else:  # n >= 2
            term1 = (4 * E + (((n - 1) ** 2) * h ** 2) / 2) * (n - 1) * rec(0, n - 1, E, h, f)
            term2 = (2 * n - 3) * rec(0, n - 2, E, h, f)
            result = (term1 - term2) / (2 * n - 1)
    elif m == 1:
        result = -1j * n * h * rec(0, n, E, h, f) / 2
    else:  # m > 1
        if m == 2:
            result = 2 * E * rec(0, n, E, h, f) - rec(0, n + 1, E, h, f) - rec(0, n - 1, E, h, f)
        else:
            term1 = (2 * E + (n * h) ** 2) * rec(m - 2, n, E, h, f)
            term2 = rec(m - 2, n + 1, E, h, f)
            term3 = rec(m - 2, n - 1, E, h, f)
            term4 = 2 * 1j * n * h * rec(m - 1, n, E, h, f)
            result = term1 - term2 - term3 - term4

    memo[(m, n, E, h, f)] = result
    return result

# Function to compute binomial-related values
def bino(m, n, p, q, E, h, f):
    result = 0
    for i in range(p + 1):
        C = math.comb(p, i)
        binomial_val = C * rec(i + m, q + n, E, h, f) * ((1j * (n + q) * h) ** (p - i))
        result += binomial_val
    return result

# Function to create the matrix for given E and f
def create_matrix(E, f, h, m, n, p, q):
    matrix = np.zeros(((m+1)*(2*n+1), (p+1)*(2*q+1)), dtype=complex)
    
    for i in range(m + 1):
        for j in range(-n, n + 1):
            for k in range(p + 1):
                for l in range(-q, q + 1):
                    matrix[i + j + n, k + l + q] = bino(i, j, k, l, E, h, f)
    
    return matrix

# Function to check if a matrix is semi-positive definite
def is_semi_pos_def(matrix):
    try:
        eigenvalues = np.linalg.eigvals(np.array(matrix))
        return np.all(eigenvalues >= 0)
    except Exception as e:
        print(f"Error in eigenvalue computation: {e}")
        return False

# Function to process a single pair of (E, f1)
def process_pair(params):
    E, f1, h, m, n, p, q = params
    memo.clear()
    A = create_matrix(E, f1, h, m, n, p, q)
    return E, f1, is_semi_pos_def(A)

def main():
    excluded_points = set()

    while True:
        # Inputs for range of E and guess value of f_1
        d = float(input("Enter range of E: "))
        f = float(input("Enter guess value of f_1: "))
        h = float(input("Enter the value of h: "))
        m = int(input("Enter range of m: "))  # Convert to int
        n = int(input("Enter range of n: "))  # Convert to int
        p = int(input("Enter range of p: "))  # Convert to int
        q = int(input("Enter range of q: "))  # Convert to int

        # Define the ranges for E and f
        range_E = np.linspace(0, d, 1000)
        range_f = np.linspace(0, f, 1000)

        # Initialize lists to store the results
        E_Values = []
        f_values = []

        # Create list of parameters to process
        params_list = [(E, f1, h, m, n, p, q) for E in range_E for f1 in range_f if (E, f1) not in excluded_points]

        # Use multiprocessing to process the pairs in parallel
        with Pool() as pool:
            results = pool.map(process_pair, params_list)

        # Collect the results
        for E, f1, is_semi_pos_def in results:
            if is_semi_pos_def:
                E_Values.append(E)
                f_values.append(f1)
            else:
                excluded_points.add((E, f1))

        # Plot the results
        plt.scatter(f_values, E_Values, cmap='viridis', edgecolor='k')
        plt.xlabel('f Values')
        plt.ylabel('E Values')
        plt.title('Scatter plot of E vs f')
        plt.colorbar(label='h Values')
        plt.show()

        # Ask the user if they want to run the program again
        cont = input("Do you want to run the program again? (yes/no): ").strip().lower()
        if cont != 'yes':
            break

if __name__ == "__main__":
    main()

