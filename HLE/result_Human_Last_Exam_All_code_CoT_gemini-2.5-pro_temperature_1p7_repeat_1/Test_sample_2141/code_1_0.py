import numpy as np
from scipy.special import genlaguerre, factorial

def solve_task():
    """
    This function calculates the ratio of the quantum to classical radial distribution
    functions for a filled n-shell of a hydrogen atom, evaluated at the radius r*
    that maximizes the classical distribution.
    
    It iterates through principal quantum numbers n = 1 to 15 to demonstrate
    the behavior of this ratio.
    """
    
    # Use atomic units (a0=1) and Z=1 for hydrogen.
    Z = 1.0
    a0 = 1.0

    def R_nl_squared(n, l, r):
        """
        Computes the square of the radial wavefunction R_nl(r).
        """
        rho = 2.0 * Z * r / (n * a0)
        
        # Normalization constant squared, from standard textbooks.
        # Use exact=True with factorial for large n to handle large integers before float conversion.
        try:
            const_sq = (2.0 * Z / (n * a0))**3 * factorial(n - l - 1, exact=True) / (2.0 * n * factorial(n + l, exact=True))
        except (ValueError, OverflowError):
            # Fallback to floating point arithmetic for very large n if needed, using log-gamma for stability.
            log_const_sq = 3*np.log(2.0/n) + gammaln(n-l) - np.log(2*n) - gammaln(n+l+1)
            const_sq = np.exp(log_const_sq)

        # Generalized Laguerre polynomial L_{n-l-1}^{2l+1}
        # scipy.special.genlaguerre(k, alpha) is L_k^alpha
        L_poly = genlaguerre(n - l - 1, 2 * l + 1)
        
        val = const_sq * np.exp(-rho) * (rho**(2 * l)) * (L_poly(rho)**2)
        return val

    def D_n(n, r):
        """
        Quantum radial distribution function for a filled n-shell.
        """
        total = 0.0
        for l in range(n):
            total += (2 * l + 1) * R_nl_squared(n, l, r)
        return r**2 * total

    def D_n_c_at_r_star(n):
        """
        Value of the classical radial distribution function at its maximum r*.
        """
        return 135.0 * Z / (128.0 * n**2)

    def ratio(n):
        """
        Computes the ratio D_n(r*)/D_n^c(r*)
        """
        # r* is the radius where D_n^c(r) is maximum
        r_s = 3.0 * n**2 / (2.0 * Z)
        
        quantum_val = D_n(n, r_s)
        classical_val = D_n_c_at_r_star(n)
        
        print(f"For n = {n}:")
        print(f"  r* = {r_s}")
        print(f"  D_n(r*) = {quantum_val}")
        print(f"  D_n_c(r*) = {classical_val}")
        print(f"  Ratio D_n(r*)/D_n_c(r*) = {quantum_val / classical_val}\n")

    # Calculate the ratio for n from 1 to 15
    for n_val in range(1, 16):
        try:
            ratio(n_val)
        except (OverflowError, ValueError) as e:
            print(f"Calculation failed for n={n_val} due to numerical instability: {e}\n")
            break
            
    print("The ratio appears to grow without bound as n increases.")
    print("Therefore, a maximum value for all positive integers n does not exist.")


solve_task()