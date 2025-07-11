import numpy as np
from scipy.special import genlaguerre, factorial

def R_nl(n, l, r, Z=1):
    """
    Calculates the normalized radial wave function R_nl(r) for a hydrogenic atom
    in atomic units (a_0 = 1).
    
    Args:
        n (int): Principal quantum number.
        l (int): Azimuthal quantum number.
        r (float): Radial distance.
        Z (int): Atomic number.
        
    Returns:
        float: Value of the radial wavefunction R_nl(r).
    """
    # rho is a dimensionless radius variable
    rho = 2 * Z * r / n
    
    # Normalization constant
    try:
        norm_factor = np.sqrt(
            (2 * Z / n)**3 * factorial(n - l - 1) / (2 * n * factorial(n + l))
        )
    except ValueError:
        # Handles cases where factorial is of a negative number, which shouldn't happen
        # for valid n, l combinations.
        return 0

    # Associated Laguerre polynomial L_{n-l-1}^{2l+1}
    laguerre_poly = genlaguerre(n - l - 1, 2 * l + 1)
    
    # Radial wave function
    return norm_factor * np.exp(-rho / 2) * (rho**l) * laguerre_poly(rho)

def D_n_quantum(n, r):
    """
    Calculates the quantum radial distribution function D_n(r).
    
    Args:
        n (int): Principal quantum number.
        r (float): Radial distance.
        
    Returns:
        float: Value of the quantum distribution D_n(r).
    """
    if r < 0:
        return 0
    
    sum_val = 0
    for l in range(n):
        R_nl_val = R_nl(n, l, r)
        sum_val += (2 * l + 1) * R_nl_val**2
        
    return r**2 * sum_val

def D_n_classical_at_r_star(n):
    """
    Calculates the value of the classical radial distribution function D_n^c(r)
    at r = r_n^* = 1.5 * n^2, where it is maximum.
    
    Args:
        n (int): Principal quantum number.
        
    Returns:
        float: Value of D_n^c(r_n^*).
    """
    return 9 / (2 * np.pi * np.sqrt(3) * n**2)

def find_max_ratio():
    """
    Finds the maximum ratio D_n(r*)/D_n^c(r*) by checking n from 1 to 30.
    """
    max_ratio = 0
    n_at_max = 0
    # Search up to n=30, which is sufficient to see the trend and find the maximum
    for n in range(1, 31):
        r_star = 1.5 * n**2
        
        d_quantum = D_n_quantum(n, r_star)
        d_classical = D_n_classical_at_r_star(n)
        
        if d_classical == 0:
            continue
        
        current_ratio = d_quantum / d_classical
        
        if current_ratio > max_ratio:
            max_ratio = current_ratio
            n_at_max = n
            
    return n_at_max, max_ratio

def main():
    """
    Main function to execute the plan and print the final answer.
    """
    # Find the n that maximizes the ratio
    n_max, _ = find_max_ratio()
    print(f"The maximum value of the ratio is found to occur at n = {n_max}.\n")
    
    # Calculate the values for n_max
    r_star_max = 1.5 * n_max**2
    d_q_max = D_n_quantum(n_max, r_star_max)
    d_c_max = D_n_classical_at_r_star(n_max)
    final_ratio = d_q_max / d_c_max

    # Print the breakdown of the calculation for n_max as per the problem's hint
    print(f"Detailed calculation for n = {n_max}:")
    print("-" * 30)

    # For n=2, r* = 1.5 * 2^2 = 6
    print(f"The classical distribution D_n^c(r) is maximized at r* = (3/2)n^2.")
    print(f"For n = {n_max}, r* = {r_star_max:.1f}\n")
    
    # For n=2, D_2(6) = 234 * exp(-6)
    num_coeff = 234 if n_max == 2 else f"D_q({r_star_max:.1f})"
    print(f"The quantum distribution at r* is D_{n_max}(r*):")
    print(f"  Analytic form for n=2: 234 * exp(-6)")
    print(f"  Numerical value: {d_q_max}\n")

    # For n=2, D_2^c(6) = 9 / (8 * pi * sqrt(3))
    den_coeff_top = 9
    den_coeff_bot = 8 if n_max == 2 else f"2*n^2"
    print(f"The classical distribution at r* is D_{n_max}^c(r*):")
    print(f"  Analytic form for n=2: {den_coeff_top} / ({den_coeff_bot} * pi * sqrt(3))")
    print(f"  Numerical value: {d_c_max}\n")

    # For n=2, Ratio = 208 * pi * sqrt(3) * exp(-6)
    ratio_coeff = 208 if n_max == 2 else "..."
    print(f"The ratio D_{n_max}(r*) / D_{n_max}^c(r*):")
    print(f"  Analytic form for n=2: {ratio_coeff} * pi * sqrt(3) * exp(-6)")
    print(f"  Maximum value: {final_ratio}\n")

    print("The maximum value of the ratio is the value calculated above.")
    print(f"{final_ratio}")


if __name__ == "__main__":
    main()