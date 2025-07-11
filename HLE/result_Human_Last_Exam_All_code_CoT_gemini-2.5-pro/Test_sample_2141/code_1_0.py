import numpy as np
from scipy.special import assoc_laguerre, factorial
import sys

def calculate_ratio(n):
    """
    Calculates the ratio D_n(r*) / D_n^c(r*) for a given integer n.
    """
    # r* is the radius that maximizes the classical distribution D_n^c(r)
    # Using np.longdouble for higher precision in intermediate calculations.
    r_star = np.longdouble(1.5 * n**2)
    
    # Calculate the normalized classical distribution value at r*.
    # D_n^c(r*) = 3*sqrt(3) / (2*pi*n^2)
    d_c_val = (3 * np.sqrt(np.longdouble(3))) / (2 * np.pi * np.longdouble(n**2))

    # --- Calculate the quantum radial distribution D_n(r*) ---
    
    # Helper function to compute the radial wavefunction R_nl(r) in atomic units
    def R_nl(n_in, l_in, r_in):
        n_ld = np.longdouble(n_in)
        l_ld = np.longdouble(l_in)
        r_ld = np.longdouble(r_in)

        # rho = 2*Z*r / (n*a0), with Z=1, a0=1
        rho = 2.0 * r_ld / n_ld
        
        # Normalization constant squared: N^2 = (2/n)^3 * (n-l-1)! / (2n * (n+l)!)
        try:
            # Use exact=False for float output, needed for large numbers
            fact_num = np.longdouble(factorial(n_in - l_in - 1, exact=False))
            fact_den = np.longdouble(factorial(n_in + l_in, exact=False))
        except ValueError: # Occurs if n-l-1 < 0
            return np.longdouble(0)
        
        const_sq = (2.0/n_ld)**3 * fact_num / (2.0 * n_ld * fact_den)
        const = np.sqrt(const_sq)
        
        # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
        # We need to compute this for potentially large arguments, which can be unstable.
        # Scipy's function should handle this reasonably for moderate n.
        lag_poly = np.longdouble(assoc_laguerre(float(rho), n_in - l_in - 1, 2 * l_in + 1))
        
        val = const * np.exp(-rho / 2.0) * rho**l_ld * lag_poly
        return val

    # D_n(r) = r^2 * sum_{l=0}^{n-1} (2l+1) * R_nl(r)^2
    total_sum = np.longdouble(0)
    for l in range(n):
        R_val = R_nl(n, l, r_star)
        total_sum += (2 * l + 1) * R_val**2
        
    d_q_val = r_star**2 * total_sum
    
    # The final ratio
    if d_c_val > 0:
        ratio = d_q_val / d_c_val
    else:
        ratio = 0.0 # Should not happen for n>0
        
    return ratio, d_q_val, d_c_val

def solve_task():
    """
    Solves the task by finding the maximum ratio over a range of n.
    """
    max_ratio = 0.0
    n_for_max = 0
    
    # The ratio appears to diverge as n increases. We will search for a maximum
    # in a reasonable range (e.g., 1 to 25).
    n_range = range(1, 26)
    
    print("This script calculates the ratio D_n(r*)/D_n^c(r*) for a given n.")
    print("r* = (3/2)n^2 is the radius that maximizes the classical probability density.")
    print("-" * 50)
    
    for n in n_range:
        try:
            current_ratio, d_q, d_c = calculate_ratio(n)
            
            # Print the components of the ratio calculation for n=2 as an example
            if n == 2:
                print("Example calculation for n=2:")
                r_star_val = 1.5 * n**2
                print(f"  r* = {r_star_val}")
                print(f"  D_2(r*) = r*^2 * [R_20(r*)^2 + 3*R_21(r*)^2] = {float(d_q):.6f}")
                print(f"  D_2^c(r*) = 3*sqrt(3) / (2*pi*n^2) = {float(d_c):.6f}")
                print(f"  Ratio for n=2: {float(d_q):.6f} / {float(d_c):.6f} = {float(current_ratio):.6f}")
                print("-" * 50)
                print("Maximum Ratio Search:")
                print("  n | Ratio D_n(r*)/D_n^c(r*)")
                print("----|-------------------------")

            if not np.isfinite(current_ratio):
                print(f" {n:2d} | Ratio became non-finite (overflow). Stopping search.")
                break

            print(f" {n:2d} | {float(current_ratio):.6f}")
            
            if current_ratio > max_ratio:
                max_ratio = current_ratio
                n_for_max = n
        except Exception as e:
            print(f"An error occurred at n={n}: {e}")
            break

    print("-" * 50)
    print("Conclusion:")
    print("The ratio appears to grow indefinitely with n, suggesting there is no finite maximum value.")
    print(f"The largest finite ratio found in the searched range (n=1 to {n_for_max}) is at n = {n_for_max}.")
    print(f"Maximum ratio found: {float(max_ratio)}")

if __name__ == '__main__':
    solve_task()