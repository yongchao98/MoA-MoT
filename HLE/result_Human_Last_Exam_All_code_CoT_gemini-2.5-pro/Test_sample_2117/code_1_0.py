import numpy as np

def solve():
    """
    Calculates the least upper bound of the product E_P*E_H * S_P*S_H.

    The derivation shows the product f(n) is:
    - f(1) = 0
    - f(n) = 2 * (1 + 1/n) for n >= 2

    We analyze this function to find its least upper bound for positive integers n.
    The sequence of values is f(1)=0, f(2)=3, f(3)=8/3, f(4)=5/2, ...
    The function is decreasing for n>=2.
    Thus, the maximum value is at n=2.
    The least upper bound of the sequence is the maximum value, which is 3.
    """
    
    n_max_val = 2
    
    # For n=2, the Cayley-Menger matrix is 3x3.
    # Eigenvalues of H are {2, -1, -1}.
    # Eigenvalues of P are {1, 1, -1}.
    
    k = n_max_val + 1
    
    # Calculate values for n=2
    E_H = (n_max_val - (-1)) / (k - 1)
    S_H = (n_max_val**2 + n_max_val * (-1)**2) / k
    
    E_P = (1 - (-1)) / (k - 1)
    S_P = 1.0 # P is orthogonal

    product = E_P * E_H * S_P * S_H

    print("Analysis for n = 2, where the maximum product occurs:")
    print(f"E_H = (2 - (-1)) / (3 - 1) = {E_H}")
    print(f"S_H = (2^2 + 2*(-1)^2) / 3 = {S_H}")
    print(f"E_P = (1 - (-1)) / (3 - 1) = {E_P}")
    print(f"S_P = 1")
    print(f"Total product = E_P * E_H * S_P * S_H = {E_P} * {E_H} * {S_P} * {S_H} = {product}")
    
    print("\nThe function f(n) = 2*(1 + 1/n) for n>=2 is a decreasing sequence.")
    print("The maximum value is at n=2, which is 3.")
    print("The least upper bound over all positive integers n is therefore 3.")

solve()
