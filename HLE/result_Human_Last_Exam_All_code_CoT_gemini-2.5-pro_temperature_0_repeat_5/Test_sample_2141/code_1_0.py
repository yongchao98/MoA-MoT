import numpy as np
from scipy.special import eval_genlaguerre
from math import sqrt, factorial, pi, exp, lgamma, log

def solve():
    """
    This script calculates the maximum value of the ratio of the quantum to classical
    radial distribution functions for a hydrogen atom's filled n-shell, evaluated
    at the radius r* that maximizes the classical distribution.
    
    The ratio is defined as: Ratio(n) = P_n(r*) / D_n^c(r*),
    where P_n(r) = D_n(r) / n^2 is the normalized quantum probability distribution,
    and D_n^c(r) is the classical probability distribution.
    """

    # Use atomic units (Bohr radius a0=1, elementary charge e=1, etc.)
    Z = 1  # Atomic number for hydrogen

    def R_nl(n, l, r):
        """
        Calculates the radial wavefunction R_nl(r) for a hydrogenic atom.
        Formula from standard quantum mechanics textbooks (e.g., Griffiths).
        R_nl(r) = N * exp(-rho/2) * rho^l * L_{n-l-1}^{2l+1}(rho)
        where rho = 2*Z*r / n
        """
        if n <= l:
            return 0
            
        rho = 2.0 * Z * r / n
        
        # Use log-gamma for factorials to handle larger n without overflow
        # log(N^2) = 3*log(2Z/n) + lgamma(n-l) - log(2n) - lgamma(n+l+1)
        try:
            # Direct computation for smaller n
            norm_factor_sq = ((2.0 * Z / n)**3 * factorial(n - l - 1)) / (2.0 * n * factorial(n + l))
        except (ValueError, OverflowError):
            # Use log-gamma for larger n
            log_norm_factor_sq = (3 * log(2.0 * Z / n) + lgamma(n - l) -
                                  log(2.0 * n) - lgamma(n + l + 1))
            norm_factor_sq = exp(log_norm_factor_sq)

        norm_factor = sqrt(norm_factor_sq)
        
        laguerre_poly = eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
        
        # The (-1) factor is a common convention, but squares away.
        return norm_factor * exp(-rho / 2.0) * (rho**l) * laguerre_poly

    def D_n_quantum(n, r):
        """
        Calculates the quantum radial distribution function D_n(r) for a filled shell.
        D_n(r) = r^2 * sum_{l=0}^{n-1} (2l+1) * |R_nl(r)|^2
        This function is normalized to n^2.
        """
        total_sum = 0
        for l in range(n):
            R_val = R_nl(n, l, r)
            total_sum += (2 * l + 1) * (R_val**2)
        return r**2 * total_sum

    def D_n_classical_peak(n):
        """
        Calculates the value of the classical radial distribution D_n^c(r) at its peak r*.
        D_n^c(r) is normalized to 1.
        The peak occurs at r* = (3/2)n^2.
        The value at the peak is D_n^c(r*) = (3*sqrt(3)) / (2*pi*n^2).
        """
        return (3 * sqrt(3)) / (2 * pi * n**2)

    max_ratio = 0
    n_at_max = 0
    details_at_max = {}

    # Loop through n to find the maximum ratio. n=30 is sufficient to show the peak.
    for n in range(1, 31):
        # r* is the radius that maximizes the classical distribution
        r_star = 1.5 * n**2
        
        # The quantum distribution D_n integrates to n^2 (the number of electrons in the shell).
        # To compare with the classical distribution D_n^c (normalized to 1), we use P_n = D_n / n^2.
        # The ratio is P_n(r*) / D_n^c(r*) = (D_n(r*) / n^2) / D_n^c(r*).
        
        d_q_val = D_n_quantum(n, r_star)
        d_c_peak_val = D_n_classical_peak(n)
        
        current_ratio = d_q_val / (n**2 * d_c_peak_val)
        
        if current_ratio > max_ratio:
            max_ratio = current_ratio
            n_at_max = n
            details_at_max = {
                "n": n,
                "r*": r_star,
                "D_n(r*)": d_q_val,
                "n^2": n**2,
                "D_n^c(r*)": d_c_peak_val,
                "Max Ratio": max_ratio
            }

    print("The maximum ratio is found for the following parameters:")
    n_val = details_at_max['n']
    r_star_val = details_at_max['r*']
    D_n_val = details_at_max['D_n(r*)']
    n_sq_val = details_at_max['n^2']
    D_nc_val = details_at_max['D_n^c(r*)']
    ratio_val = details_at_max['Max Ratio']
    
    print(f"The ratio is D_n(r*) / (n^2 * D_n^c(r*))")
    print(f"Principal quantum number n = {n_val}")
    print(f"Classical peak radius r* = {r_star_val:.4f}")
    print(f"Quantum distribution at r*, D_n(r*) = {D_n_val:.4e}")
    print(f"Normalization factor, n^2 = {n_sq_val}")
    print(f"Classical distribution at r*, D_n^c(r*) = {D_nc_val:.4e}")
    print(f"Maximum ratio value = {ratio_val:.4f}")

solve()