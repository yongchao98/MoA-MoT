import numpy as np
from scipy.special import genlaguerre
from math import factorial

def R_nl(n, l, r):
    """
    Calculates the value of the hydrogenic radial wavefunction R_nl(r)
    in atomic units (Z=1, a_0=1).
    """
    # Normalization constant for R_nl
    try:
        norm_const = np.sqrt((2.0/n)**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    except ValueError:
        # This can happen if n-l-1 < 0, but our loop ensures l < n.
        return 0.0

    # Argument for the functions
    rho = 2.0 * r / n
    
    # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    
    # Full radial wavefunction
    wavefunction_val = norm_const * np.exp(-rho / 2.0) * (rho**l) * laguerre_poly(rho)
    return wavefunction_val

def D_n(n, r):
    """
    Calculates the quantum radial distribution function D_n(r) for a filled n-shell.
    """
    d_sum = 0.0
    for l in range(n):
        r_nl_val = R_nl(n, l, r)
        d_sum += (2 * l + 1) * r_nl_val**2
            
    return r**2 * d_sum

def D_nc_max(n):
    """
    Calculates the maximum value of the classical radial distribution function, D_n^c(r^*).
    """
    return (3 * np.sqrt(3)) / (np.pi * n**2)

def solve_and_print():
    """
    Finds the maximum ratio D_n(r*)/D_n^c(r*) by checking integer values of n,
    and prints the detailed calculation for the n that gives the maximum.
    """
    max_ratio = 0.0
    n_for_max_ratio = 0

    # The ratio approaches 1 for large n, so we only need to check small n.
    for n_test in range(1, 16):
        r_star_test = 1.5 * n_test**2
        d_n_val_test = D_n(n_test, r_star_test)
        d_nc_val_test = D_nc_max(n_test)
        ratio = d_n_val_test / d_nc_val_test
        if ratio > max_ratio:
            max_ratio = ratio
            n_for_max_ratio = n_test

    # Now, perform and print the calculation for the n that yielded the max ratio.
    n = n_for_max_ratio
    r_star = 1.5 * n**2
    
    d_n_val_at_max = D_n(n, r_star)
    d_nc_val_at_max = D_nc_max(n)
    final_ratio = d_n_val_at_max / d_nc_val_at_max
    
    print(f"The analysis shows the maximum ratio occurs for n = {n}.")
    print("\n--- Calculation for n = {} ---".format(n))
    print("The radius r* that maximizes the classical distribution D_n^c(r) is:")
    print(f"r* = (3/2) * n^2 = (3/2) * {n}^2 = {r_star}")
    
    print("\nThe maximum value of the classical distribution D_n^c(r) is:")
    print(f"D_{n}^c(r*) = (3 * sqrt(3)) / (pi * n^2) = (3 * sqrt(3)) / (pi * {n}^2) = {d_nc_val_at_max}")

    print("\nThe value of the quantum distribution D_n(r) at r* is:")
    print(f"D_{n}(r*) = (r*)^2 * Sum[(2l+1)*|R_nl(r*)|^2, l=0..n-1] = {d_n_val_at_max}")

    print("\nThe ratio D_n(r*)/D_n^c(r*) is therefore:")
    print(f"{d_n_val_at_max} / {d_nc_val_at_max} = {final_ratio}")

solve_and_print()