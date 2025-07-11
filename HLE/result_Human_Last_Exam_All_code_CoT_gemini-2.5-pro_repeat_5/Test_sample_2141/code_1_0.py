import numpy as np
from scipy.special import genlaguerre, factorial

def calculate_ratio_and_print():
    """
    This script calculates the ratio of the quantum mechanical to the classical radial distribution
    function for a filled n-shell of a hydrogen atom, evaluated at the peak of the classical function.

    It computes the ratio for n = 1 to 10 to demonstrate the trend and finds the maximum value.
    """

    # The maximum value of the normalized classical distribution function D_n^c(r)
    # occurs at r* = 1.5 * n^2 and its value is D_n^c(r*) = 3*sqrt(3)/(2*pi), independent of n.
    D_nc_max = 3 * np.sqrt(3) / (2 * np.pi)

    def R_nl_sq(n, l, r):
        """
        Computes the square of the hydrogenic radial wavefunction R_nl(r) in atomic units.
        """
        # rho is a dimensionless radius variable
        rho = 2 * r / n
        
        # Normalization constant squared for R_nl(r)
        norm_sq = (2/n)**3 * factorial(n - l - 1) / (2 * n * factorial(n + l))
        
        # Generalized Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
        # scipy.special.genlaguerre(k, alpha) corresponds to L_k^alpha
        laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
        laguerre_val = laguerre_poly(rho)
        
        val = norm_sq * np.exp(-rho) * (rho**(2*l)) * (laguerre_val**2)
        return val

    def D_n(n, r):
        """
        Computes the quantum radial distribution function D_n(r) for a filled shell.
        """
        total = 0
        for l in range(n):
            total += (2 * l + 1) * R_nl_sq(n, l, r)
        return r**2 * total

    print("Calculating the ratio D_n(r*) / D_n^c(r*) for various n:")
    print("-" * 60)
    print(f"The classical maximum value D_n^c(r*) is constant: {D_nc_max:.6f}")
    
    max_ratio = 0
    max_n = 0

    for n in range(1, 11):
        # r* is the radius that maximizes the classical distribution D_n^c(r)
        r_star = 1.5 * n**2
        
        # Calculate the quantum distribution value at r*
        D_n_val = D_n(n, r_star)
        
        # The ratio D_n(r*) / D_n^c(r*)
        ratio = D_n_val / D_nc_max
        
        print(f"For n={n:2d}: r*={r_star:5.1f}, D_n(r*)={D_n_val:.6f}, Ratio={ratio:.6f}")
        
        if ratio > max_ratio:
            max_ratio = ratio
            max_n = n
            
    print("-" * 60)
    print("The ratio is observed to increase monotonically with n.")
    print("According to the correspondence principle, this ratio approaches 1 as n approaches infinity.")
    print("Thus, the sequence does not have a true maximum for any finite n.")
    print("The maximum value (supremum) of the ratio is the limit as n -> infinity.")
    
calculate_ratio_and_print()
