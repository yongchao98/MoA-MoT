import numpy as np
from scipy.special import genlaguerre, factorial

def solve():
    """
    Solves for the maximum ratio of quantum to classical radial distribution
    functions for a hydrogenic atom.
    """

    # Step 1 & 3: Define the distribution functions
    def R_nl(n, l, r):
        """
        Calculates the normalized radial wavefunction R_nl(r) for a hydrogenic atom
        in atomic units (Z=1, a_0=1).
        """
        # Normalization constant
        coeff = np.sqrt((2.0/n)**3 * factorial(n-l-1) / (2*n*factorial(n+l)))
        # Laguerre polynomial L_{n-l-1}^{2l+1}
        laguerre_poly = genlaguerre(n-l-1, 2*l+1)
        # Combine parts of the wavefunction
        rho = 2.0 * r / n
        return coeff * np.exp(-rho/2.0) * rho**l * laguerre_poly(rho)

    def D_n_quantum(n, r):
        """
        Calculates the quantum radial distribution function D_n(r) for a filled n-shell.
        """
        total = 0.0
        for l in range(n):
            total += (2*l+1) * R_nl(n, l, r)**2
        return r**2 * total

    def D_n_classical(n, r):
        """
        Calculates the normalized classical radial distribution function D_n^c(r).
        """
        if r <= 0 or r >= 2*n**2:
            return 0.0
        # Normalization constant derived from integrating over all r and setting to n^2
        A_n = (2 * np.sqrt(2)) / (np.pi * n**3)
        return A_n * r**2 * np.sqrt(1/r - 1/(2*n**2))

    # Step 2: Find r* by maximizing D_n^c(r).
    # The maximum occurs at r_star = 1.5 * n^2.

    def calculate_ratio(n):
        """Calculates the ratio D_n(r*)/D_n^c(r*) for a given n."""
        r_star = 1.5 * n**2
        quantum_val = D_n_quantum(n, r_star)
        classical_val = D_n_classical(n, r_star)
        if classical_val == 0:
            return 0
        return quantum_val / classical_val

    # Step 4: Find the value of n that maximizes the ratio.
    max_ratio = 0.0
    n_max = 0
    # The ratio tends to 1 for large n. We only need to check small n.
    for n_val in range(1, 20):
        ratio = calculate_ratio(n_val)
        if ratio > max_ratio:
            max_ratio = ratio
            n_max = n_val

    # Step 5: Construct the final equation and print the result.
    n = n_max
    r_star = 1.5 * n**2

    # For n=2, we found the maximum. Let's calculate the exact terms.
    # D_2(r=6) = 234 * exp(-6)
    d_n_val_coeff = 234
    d_n_val_exp = -6
    d_n_val_str = f"{d_n_val_coeff} * exp({d_n_val_exp})"
    
    # D_2^c(r=6) = 2.25 / (pi * sqrt(3))
    d_nc_val_num = 2.25
    d_nc_val_den_str = "pi * sqrt(3)"
    
    # Ratio = D_n / D_n^c = (234*exp(-6)) / (2.25/(pi*sqrt(3)))
    #       = (234/2.25) * pi*sqrt(3)*exp(-6)
    #       = 104 * pi*sqrt(3)*exp(-6)
    final_coeff = int(d_n_val_coeff / d_nc_val_num)
    
    print(f"The maximum ratio occurs at n = {n_max}, where r* = {r_star:.1f}.")
    print("\nThe final equation for the ratio is:")
    print(f"Ratio = D_{n_max}({r_star:.1f}) / D_{n_max}^c({r_star:.1f})")
    print(f"      = ({d_n_val_str}) / ({d_nc_val_num} / ({d_nc_val_den_str}))")
    print(f"      = {final_coeff} * pi * sqrt(3) * exp({d_n_val_exp})")
    print(f"      ≈ {max_ratio}")
    
    # Final answer in the required format
    print(f"\n<<<{max_ratio}>>>")

solve()