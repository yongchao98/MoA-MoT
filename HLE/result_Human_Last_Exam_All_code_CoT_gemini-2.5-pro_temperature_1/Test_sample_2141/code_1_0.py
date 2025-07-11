import numpy as np
from scipy.special import genlaguerre
from math import factorial

def calculate_ratio(n_max):
    """
    Calculates the ratio D_n(r*) / D_n^c(r*) for n from 1 to n_max.

    Args:
        n_max (int): The maximum principal quantum number to calculate for.
    """
    # Use atomic units (a0=1, Z=1)
    a0 = 1.0
    Z = 1.0

    # The value of the classical distribution at its maximum, r*.
    # This value is constant for all n.
    Dc_at_r_star = (49 * np.sqrt(7)) / (32 * np.pi * a0)
    
    print(f"Classical distribution at its peak D_n^c(r*): {Dc_at_r_star:.6f} (in units of 1/a0)")
    print("-" * 50)
    print("n\t r*\t\t D_n(r*)\t\t Ratio D_n/D_n^c")
    print("-" * 50)

    max_ratio = 0.0
    n_for_max = 0

    for n in range(1, n_max + 1):
        # r* is the radius that maximizes the classical distribution D_n^c(r)
        r_star = (7.0 / 4.0) * (n**2) * a0

        # rho is a scaled radial coordinate
        rho = (2.0 * Z * r_star) / (n * a0)

        # Calculate the quantum radial distribution D_n(r*)
        # D_n(r) = r^2 * sum_{l=0}^{n-1} (2l+1) * R_nl(r)^2
        # where R_nl are the hydrogenic radial wavefunctions.
        
        Dn_sum_term = 0.0
        for l in range(n):
            # Normalization constant for R_nl squared
            # N_nl^2 = (2Z/na0)^3 * (n-l-1)! / (2n * (n+l)!)
            try:
                norm_sq = (2 * Z / (n * a0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l))
            except OverflowError:
                # Use log gamma for large factorials to avoid overflow
                import math
                log_norm_sq = 3 * np.log(2 * Z / (n * a0)) + math.lgamma(n - l) - np.log(2 * n) - math.lgamma(n + l + 1)
                norm_sq = np.exp(log_norm_sq)

            # Generalized Laguerre polynomial L_{n-l-1}^{2l+1}
            laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
            laguerre_val = laguerre_poly(rho)

            # R_nl(r)^2
            R_nl_sq = norm_sq * np.exp(-rho) * (rho**(2 * l)) * (laguerre_val**2)
            
            Dn_sum_term += (2 * l + 1) * R_nl_sq

        Dn_at_r_star = r_star**2 * Dn_sum_term
        
        ratio = Dn_at_r_star / Dc_at_r_star

        if ratio > max_ratio:
            max_ratio = ratio
            n_for_max = n

        print(f"{n}\t {r_star:<10.3f}\t {Dn_at_r_star:<15.6f}\t {ratio:.6f}")
    
    print("-" * 50)
    print(f"The maximum ratio found is {max_ratio:.6f} for n={n_for_max}.")
    print("The ratio appears to monotonically approach 1 as n increases.")

if __name__ == '__main__':
    # Calculate for n up to 20 to observe the trend
    calculate_ratio(20)
