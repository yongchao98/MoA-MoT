import numpy as np
from scipy.special import genlaguerre
from math import factorial

def get_classical_dist_max_val(n):
    """
    Calculates the value of the normalized classical radial distribution D_n^c(r)
    at the radius r* that maximizes it. Atomic units are used (a_0=1, e=1, etc.).
    The radius r* that maximizes D_n^c(r) is r* = 1.5 * n^2.
    The value of the distribution at this point is D_n^c(r*) = (3 * sqrt(3)) / (2 * pi * n^3).
    
    Args:
        n (int): The principal quantum number.
        
    Returns:
        float: The value of D_n^c(r*).
    """
    if n <= 0:
        return 0
    return (3 * np.sqrt(3)) / (2 * np.pi * n**3)

def get_radial_wavefunction(n, l, r):
    """
    Calculates the hydrogenic radial wavefunction R_nl(r) in atomic units (Z=1, a_0=1).
    
    Args:
        n (int): Principal quantum number.
        l (int): Azimuthal quantum number.
        r (float): Radial distance.
        
    Returns:
        float: The value of R_nl(r).
    """
    # Normalization constant for the radial wavefunction
    coeff = np.sqrt((2.0/n)**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    
    # The argument for the polynomials and exponential
    rho = 2.0 * r / n
    
    # Generalized Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
    # scipy.special.genlaguerre is defined as L_k^(alpha)
    k = n - l - 1
    alpha = 2 * l + 1
    laguerre_poly = genlaguerre(k, alpha)
    
    return coeff * np.exp(-rho / 2.0) * rho**l * laguerre_poly(rho)

def get_quantum_dist_val(n, r):
    """
    Calculates the value of the quantum radial distribution for a filled n-shell,
    D_n(r) = r^2 * sum_{l=0}^{n-1} (2l+1) * R_nl(r)^2.
    
    Args:
        n (int): Principal quantum number.
        r (float): Radial distance.
        
    Returns:
        float: The value of D_n(r).
    """
    total_sum = 0
    for l in range(n):
        R_nl_val = get_radial_wavefunction(n, l, r)
        total_sum += (2*l + 1) * (R_nl_val**2)
    
    return r**2 * total_sum

def find_max_ratio():
    """
    Calculates the ratio D_n(r*)/D_n^c(r*) for a range of n values,
    finds the maximum ratio, and prints the results.
    """
    max_ratio = 0
    n_max = 0
    r_star_at_max = 0
    d_n_at_max = 0
    d_nc_at_max = 0

    # Iterate through n from 1 to 30 to find the maximum
    for n in range(1, 31):
        r_star = 1.5 * n**2
        d_n_val = get_quantum_dist_val(n, r_star)
        d_nc_val = get_classical_dist_max_val(n)
        
        if d_nc_val == 0:
            continue
            
        current_ratio = d_n_val / d_nc_val
        
        if current_ratio > max_ratio:
            max_ratio = current_ratio
            n_max = n
            r_star_at_max = r_star
            d_n_at_max = d_n_val
            d_nc_at_max = d_nc_val

    print("Analysis complete. The maximum ratio is found for n =", n_max)
    print("For n = {}:".format(n_max))
    print("r* (radius of max classical probability) = {:.4f}".format(r_star_at_max))
    print("D_n(r*) (Quantum value) = {:.4f}".format(d_n_at_max))
    print("D_n^c(r*) (Classical value) = {:.4f}".format(d_nc_at_max))
    print("Maximum Ratio D_n(r*)/D_n^c(r*) = {:.4f}".format(max_ratio))

if __name__ == '__main__':
    find_max_ratio()