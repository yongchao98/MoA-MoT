import numpy as np
from scipy.special import genlaguerre
from math import factorial, sqrt, pi, exp

def R_nl_squared(n, l, r):
    """
    Calculates the square of the radial wavefunction R_nl(r) for a hydrogenic atom.
    Uses atomic units (a_0=1, Z=1).
    """
    # rho is a scaled radius variable
    rho = 2.0 * r / n
    
    # Normalization constant squared
    try:
        norm_sq = (2.0/n)**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l))
    except ValueError:
        # factorial of a negative number, occurs for l >= n, should not happen in our loop
        return 0
    
    # Generalized Laguerre polynomial L_{n-l-1}^{2l+1}(\rho)
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    laguerre_val_sq = laguerre_poly(rho)**2
    
    return norm_sq * exp(-rho) * (rho**(2*l)) * laguerre_val_sq

def D_n(n, r):
    """
    Calculates the quantum radial distribution function D_n(r) for a filled n-shell.
    This function is not normalized to 1; its integral over r is n^2.
    """
    sum_val = 0
    for l in range(n):
        sum_val += (2*l + 1) * R_nl_squared(n, l, r)
    return r**2 * sum_val

def get_ratio(n):
    """
    Calculates the ratio D_n(r*)/(n^2 * D_n^c(r*)).
    """
    # r* is the radius that maximizes the classical distribution D_n^c(r)
    r_star = 1.5 * n**2
    
    # Calculate the quantum distribution value at r*
    dn_val = D_n(n, r_star)
    
    # Normalize the quantum distribution value
    dn_val_normalized = dn_val / n**2
    
    # Calculate the classical distribution value at r*
    # The analytical expression for D_n_c(r*) is (3*sqrt(3))/(2*pi*n^2)
    dnc_val = (3 * sqrt(3)) / (2 * pi * n**2)
    
    # Return the ratio, and the components for detailed output
    return dn_val_normalized / dnc_val, dn_val_normalized, dnc_val, r_star

def solve():
    """
    Finds the maximum value of the ratio for n from 1 to 20.
    """
    max_ratio = 0
    n_max = 0
    
    # Store all results to print them later
    results = {}

    # We check up to n=20, which is sufficient to see the trend
    # as the ratio approaches 1 for large n.
    for n in range(1, 21):
        ratio_val, dn_norm, dnc, r_s = get_ratio(n)
        results[n] = {'ratio': ratio_val, 'dn_norm': dn_norm, 'dnc': dnc, 'r_star': r_s}
        if ratio_val > max_ratio:
            max_ratio = ratio_val
            n_max = n
            
    print("This script calculates the ratio of the normalized quantum radial density to the classical radial density, D_n(r*)/(n^2 * D_n^c(r*)), evaluated at the radius r* that maximizes the classical density.\n")
    print(f"The maximum value of the ratio occurs at n = {n_max}.\n")
    
    res_max = results[n_max]
    n = n_max
    dn_norm_max = res_max['dn_norm']
    dnc_max = res_max['dnc']
    r_star_max = res_max['r_star']
    
    print(f"Calculation for n = {n}:")
    print(f"The radius r* that maximizes D_n^c(r) is (3/2)*n^2 = {r_star_max:.4f}")
    
    print(f"The value of the classical density at r* is D_n^c(r*) = (3*sqrt(3))/(2*pi*n^2) = {dnc_max:.6f}")
    
    print(f"The value of the normalized quantum density at r* is D_n(r*)/n^2 = {dn_norm_max:.6f}")
    
    print("\nThe ratio is therefore:")
    print(f"{dn_norm_max:.6f} / {dnc_max:.6f} = {max_ratio:.6f}\n")
    
    print("The maximum value of the ratio is:")
    print(max_ratio)
    
solve()