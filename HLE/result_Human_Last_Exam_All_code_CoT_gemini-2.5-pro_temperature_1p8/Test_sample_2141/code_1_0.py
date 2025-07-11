import numpy as np
from scipy.special import genlaguerre, factorial

def R_nl(n, l, r):
    """
    Calculates the normalized radial wavefunction R_nl(r) for a hydrogenic atom.
    Assumes atomic units (Z=1, a_0=1).
    """
    # The variable rho in the standard formula
    rho = 2.0 * r / n
    
    # Normalization constant
    # Note: factorial from scipy.special handles arrays, np.math.factorial does not.
    norm_factor = np.sqrt(
        (2.0 / n)**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l))
    )
    
    # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(\rho)
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    
    # The wavefunction expression
    return norm_factor * np.exp(-rho / 2.0) * (rho**l) * laguerre_poly(rho)

def D_n_quantum(n, r):
    """
    Calculates the quantum radial distribution function D_n(r) for a filled shell.
    """
    total_sum = 0
    for l in range(n):
        R_val = R_nl(n, l, r)
        total_sum += (2 * l + 1) * R_val**2
    return r**2 * total_sum

def solve_problem():
    """
    Calculates the ratio D_n(r*)/D_n^c(r*) for various n and finds the maximum.
    """
    # Value of the classical distribution at its peak is constant
    D_classical_peak = (3 * np.sqrt(3)) / (2 * np.pi)

    max_ratio = 0.0
    n_at_max = 0

    print("--- Calculating the ratio of Quantum to Classical Radial Distributions ---")
    print(f"Classical peak value D_n^c(r*) = {D_classical_peak:.6f} (constant for all n)")
    print("-" * 60)
    print("n      r*         D_n(r*)        Ratio D_n(r*)/D_n^c(r*)")
    print("-" * 60)

    # Calculate for n=1 to 10, which should be sufficient to see the trend and find the max
    for n in range(1, 11):
        # r* is the radius where the classical distribution is maximal
        r_star = 1.5 * n**2
        
        # Calculate the quantum distribution at r*
        D_quantum_at_r_star = D_n_quantum(n, r_star)
        
        # Compute the ratio
        current_ratio = D_quantum_at_r_star / D_classical_peak
        
        print(f"{n:<5}  {r_star:<10.2f} {D_quantum_at_r_star:<14.6f} {current_ratio:.6f}")

        if current_ratio > max_ratio:
            max_ratio = current_ratio
            n_at_max = n
    
    print("-" * 60)
    print(f"The maximum ratio occurs at n = {n_at_max}")
    print(f"The maximum value of the ratio is {max_ratio:.6f}")

if __name__ == '__main__':
    solve_problem()
