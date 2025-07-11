import numpy as np
from scipy.special import comb, factorial

def laguerre_poly(k, alpha, x):
    """
    Calculates the generalized Laguerre polynomial L_k^(alpha)(x) using the
    summation formula, which corresponds to the definition widely used
    in mathematical physics textbooks.
    """
    if k < 0:
        return 0
    total = 0
    for i in range(k + 1):
        # Using SciPy's combination (nCk) and factorial functions for precision
        term = (-1)**i * comb(k + alpha, k - i, exact=True) * (x**i / factorial(i, exact=True))
        total += term
    return total

def R_nl(n, l, r):
    """
    Calculates the Hydrogenic radial wavefunction R_nl(r) in atomic units (Z=1, a_0=1).
    """
    # The variable 'rho' in the standard wavefunction formula
    rho = 2.0 * r / n
    
    # Prefactor for normalization
    # Using np.math.factorial for large numbers before taking sqrt
    fact_ratio = np.math.factorial(n - l - 1) / (2.0 * n * np.math.factorial(n + l))
    prefactor = np.sqrt((2.0/n)**3 * fact_ratio)
    
    # The Laguerre polynomial part, L_{n-l-1}^{2l+1}(\rho)
    lag_poly = laguerre_poly(n - l - 1, 2*l + 1, rho)
    
    return prefactor * np.exp(-rho/2.0) * (rho**l) * lag_poly

def D_n_quantum(n):
    """
    Calculates the quantum radial distribution function D_n(r) at the
    point r* where the classical distribution is maximal.
    r* = 1.5 * n^2
    """
    r_star = 1.5 * n**2
    
    total = 0.0
    # Sum over all orbital angular momenta l from 0 to n-1
    for l in range(n):
        R_nl_val = R_nl(n, l, r_star)
        total += (2*l + 1) * R_nl_val**2
            
    return r_star**2 * total

def find_max_ratio():
    """
    Calculates the ratio D_n(r*)/D_c(r*) for n=1, 2, 3, ...
    and finds the maximum value.
    """
    # The value of the normalized classical distribution at its maximum is constant
    D_n_classical_max = (3 * np.sqrt(3)) / (2 * np.pi)
    
    max_ratio = 0.0
    n_at_max = 0
    
    print("Calculating the ratio for different values of n:")
    print(" n |     Ratio")
    print("---|------------")
    # Loop over a range of n sufficient to find the peak. The ratio converges to 1.
    for n in range(1, 21):
        D_n_q_at_max = D_n_quantum(n)
        ratio = D_n_q_at_max / D_n_classical_max
        print(f"{n:2d} | {ratio:10.6f}")
        if ratio > max_ratio:
            max_ratio = ratio
            n_at_max = n
    
    # Outputting the numbers for the final equation as requested
    final_D_n_q = D_n_quantum(n_at_max)
    
    print("\n--- Details of the Maximum Ratio Calculation ---")
    print(f"The maximum ratio is found to occur for n = {n_at_max}.")
    print(f"For n = {n_at_max}, the classical distribution is maximized at r* = {1.5*n_at_max**2}.")
    
    # Final equation: Max Ratio = D_n(r*) / D_c(r*)
    print("\nFinal Equation:")
    print(f"Numerator: Quantum value D_n(r*) = {final_D_n_q}")
    print(f"Denominator: Classical value D_c(r*) = {D_n_classical_max}")
    print(f"Result: Maximum Ratio = {final_D_n_q} / {D_n_classical_max} = {max_ratio}")
    

if __name__ == "__main__":
    find_max_ratio()
