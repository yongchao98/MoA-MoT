import numpy as np
from scipy.special import genlaguerre
from math import factorial, sqrt, pi, exp

def get_R_nl_at_r(n, l, r):
    """
    Computes the value of the hydrogenic radial wavefunction R_nl(r).
    Assumes atomic units (Z=1, a0=1).
    """
    rho = 2.0 * r / n
    
    # Normalization constant from standard formula
    try:
        norm_const_sq = (2.0 / n)**3 * (factorial(n - l - 1)) / (2.0 * n * factorial(n + l))
        norm_const = sqrt(norm_const_sq)
    except (ValueError, OverflowError):
        # Use log-gamma for large n to avoid overflow
        from scipy.special import loggamma
        log_norm_sq = 3 * np.log(2.0 / n) + loggamma(n - l) - np.log(2.0 * n) - loggamma(n + l + 1)
        norm_const = np.exp(0.5 * log_norm_sq)

    # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    
    # Combine all parts of the wavefunction
    val = norm_const * np.exp(-rho / 2.0) * (rho**l) * laguerre_poly(rho)
    return val

def main():
    """
    Calculates the ratio for n=1 to 20 and finds the maximum.
    """
    max_ratio = 0.0
    n_at_max = 0
    max_dn_val = 0
    max_dnc_val = 0
    
    print("This script calculates the ratio of the quantum to classical radial probability")
    print("densities for a hydrogen atom's filled n-shell, evaluated at the peak of the")
    print("classical distribution.\n")
    print("Ratio = (D_n(r*) / n^2) / D_n^c(r*)")
    print("--------------------------------------------------")
    print(" n | r*       | D_n(r*)/n^2  | D_n^c(r*)    | Ratio")
    print("--------------------------------------------------")

    # Loop through a range of n values
    for n in range(1, 21):
        r_star = 1.5 * n**2
        rho_star = 3.0 * n

        # --- Calculate quantum distribution value D_n(r*) ---
        dn_sum = 0.0
        for l in range(n):
            R_val_sq = get_R_nl_at_r(n, l, r_star)**2
            dn_sum += (2 * l + 1) * R_val_sq
        
        dn_val = r_star**2 * dn_sum
        
        # Normalize per state
        dn_per_state = dn_val / n**2

        # --- Calculate classical distribution value D_n^c(r*) ---
        # Value at maximum is known analytically
        dnc_val = (3.0 * sqrt(3.0)) / (2.0 * pi * n**2)

        # Calculate the final ratio
        if dnc_val == 0:
            ratio = float('inf')
        else:
            ratio = dn_per_state / dnc_val
        
        print(f"{n:2d} | {r_star:<8.1f} | {dn_per_state:<12.6f} | {dnc_val:<12.6f} | {ratio:.6f}")

        # Keep track of the maximum ratio found
        if ratio > max_ratio:
            max_ratio = ratio
            n_at_max = n
            max_dn_val = dn_per_state
            max_dnc_val = dnc_val
    
    print("--------------------------------------------------")
    print("\nCalculation complete.")
    print(f"The maximum value of the ratio occurs at n = {n_at_max}.")
    print(f"At this point, D_n(r*)/n^2 is {max_dn_val:.6f}.")
    print(f"The classical value D_n^c(r*) is {max_dnc_val:.6f}.")
    print(f"The maximum ratio is D_n(r*)/n^2 / D_n^c(r*) = {max_dn_val:.6f} / {max_dnc_val:.6f} = {max_ratio:.6f}.")

if __name__ == "__main__":
    main()