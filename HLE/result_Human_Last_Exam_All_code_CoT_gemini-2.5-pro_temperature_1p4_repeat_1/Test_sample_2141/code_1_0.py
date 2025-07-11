import math
from scipy.special import genlaguerre

def get_R_nl_squared(n, l, r):
    """
    Calculates the square of the radial wavefunction |R_nl(r)|^2 for a hydrogenic atom.
    Formula is in atomic units (Z=1, a0=1).
    """
    # rho is a dimensionless radius variable
    rho = 2 * r / n
    
    # Prefactor term involving factorials, computed using lgamma for numerical stability
    log_prefactor_A = 3 * (math.log(2) - math.log(n))
    log_prefactor_B = math.lgamma(n - l) - (math.log(2 * n) + math.lgamma(n + l + 1))
    prefactor = math.exp(log_prefactor_A + log_prefactor_B)
    
    # Laguerre polynomial term
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    laguerre_val = laguerre_poly(rho)
    
    # Combine all parts to get |R_nl(r)|^2
    r_nl_sq = prefactor * (rho**(2 * l)) * math.exp(-rho) * (laguerre_val**2)
    return r_nl_sq

def get_D_n(n, r):
    """
    Calculates the quantum radial distribution function D_n(r) for a filled n-shell.
    """
    d_n = 0
    for l in range(n):
        r_nl_sq = get_R_nl_squared(n, l, r)
        d_n += (2 * l + 1) * r_nl_sq
    return r**2 * d_n

def solve():
    """
    Finds the maximum value of the ratio D_n(r*)/D_n^c(r*) and prints the result.
    """
    print("This script finds the maximum value of the ratio D_n(r*)/D_n^c(r*).")
    print("We will evaluate the ratio for n = 1, 2, 3, ... to find the maximum.\n")
    
    max_ratio = 0
    max_n = 0
    
    # Check n from 1 to 10, which is sufficient to find the maximum
    for n_val in range(1, 11):
        # r* is the radius that maximizes the classical distribution D_n^c(r)
        r_star = 1.5 * n_val**2
        
        # Calculate the quantum distribution at r*
        d_n_at_r_star = get_D_n(n_val, r_star)
        
        # The maximum value of the classical distribution D_n^c(r*) is a constant
        d_nc_at_r_star = (3 * math.sqrt(3)) / (2 * math.pi)
        
        # Calculate the ratio
        ratio = d_n_at_r_star / d_nc_at_r_star
        
        # print(f"For n = {n_val}, Ratio = {ratio:.6f}")
        
        if ratio > max_ratio:
            max_ratio = ratio
            max_n = n_val

    print(f"The maximum value of the ratio is found at n = {max_n}.")
    print("-" * 40)
    print(f"Calculation for n = {max_n}:")

    # Recalculate for the max_n to show details
    r_star = 1.5 * max_n**2
    d_n_at_r_star = get_D_n(max_n, r_star)
    d_nc_at_r_star = (3 * math.sqrt(3)) / (2 * math.pi)

    print(f"The classical distribution D_n^c(r) is maximized at r* = (3/2)n^2.")
    print(f"For n = {max_n}, r* = {r_star}")
    
    print(f"The value of the quantum distribution at r* is D_{max_n}(r*) = {d_n_at_r_star:.6f}")
    print(f"The value of the classical distribution at r* is D_{max_n}^c(r*) = 3*sqrt(3)/(2*pi) = {d_nc_at_r_star:.6f}")
    
    print("\nThe ratio is D_n(r*)/D_n^c(r*).")
    print(f"Maximum ratio = {d_n_at_r_star:.6f} / {d_nc_at_r_star:.6f}")
    
    print(f"Final Answer = {max_ratio}")
    print("The exact analytical expression for this value is 52 * sqrt(3) * pi * e^(-6).")

solve()
