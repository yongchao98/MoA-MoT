import numpy as np

def get_R_nl_sq(n, l):
    """
    Returns a function for the square of the radial wavefunction R_nl(r)
    for given n and l, using formulas from standard texts (e.g., Wikipedia).
    These formulas are for a hydrogenic atom with Z=1 in atomic units (a_0=1).
    """
    if n == 1 and l == 0:
        return lambda r: (2 * np.exp(-r))**2
    elif n == 2 and l == 0:
        return lambda r: (1/(2*np.sqrt(2)) * (2-r) * np.exp(-r/2))**2
    elif n == 2 and l == 1:
        return lambda r: (1/(2*np.sqrt(6)) * r * np.exp(-r/2))**2
    elif n == 3 and l == 0:
        return lambda r: (2/(81*np.sqrt(3)) * (27 - 18*r + 2*r**2) * np.exp(-r/3))**2
    elif n == 3 and l == 1:
        # This is the correctly normalized version
        return lambda r: ((4*np.sqrt(2))/(9*np.sqrt(3)) * r * (1 - r/6) * np.exp(-r/3))**2
    elif n == 3 and l == 2:
        # This is the correctly normalized version
        return lambda r: ((4/(81*np.sqrt(30))) * r**2 * np.exp(-r/3))**2
    else:
        # Return a function that is always zero for unsupported n, l
        return lambda r: 0

def calculate_quantum_distribution_d_n(n, r):
    """
    Calculates the single-particle quantum radial distribution d_n(r) = D_n(r)/n^2.
    """
    D_n_r = 0
    for l in range(n):
        R_nl_sq_func = get_R_nl_sq(n, l)
        D_n_r += (2*l + 1) * R_nl_sq_func(r)
    D_n_r *= r**2
    return D_n_r / n**2

def calculate_classical_distribution_at_peak(n):
    """
    Calculates the value of the classical distribution D_n^c(r) at its peak r* = 1.5*n^2.
    """
    return (3 * np.sqrt(3)) / (2 * np.pi * n**2)

def solve():
    """
    Calculates the ratio for n=1, 2, 3 and finds the maximum.
    """
    max_ratio = 0
    max_n = 0
    
    for n in range(1, 4):
        r_star = 1.5 * n**2
        
        d_n_at_r_star = calculate_quantum_distribution_d_n(n, r_star)
        D_n_c_at_r_star = calculate_classical_distribution_at_peak(n)
        
        ratio = d_n_at_r_star / D_n_c_at_r_star
        
        if ratio > max_ratio:
            max_ratio = ratio
            max_n = n
            
    print(f"The maximum value of the ratio is {max_ratio:.4f}, which occurs for n = {max_n}.")
    print(f"The calculation for n=1 gives a ratio of {calculate_quantum_distribution_d_n(1, 1.5*1**2)/calculate_classical_distribution_at_peak(1):.4f}")
    print(f"The calculation for n=2 gives a ratio of {calculate_quantum_distribution_d_n(2, 1.5*2**2)/calculate_classical_distribution_at_peak(2):.4f}")
    print(f"The calculation for n=3 gives a ratio of {calculate_quantum_distribution_d_n(3, 1.5*3**2)/calculate_classical_distribution_at_peak(3):.4f}")
    
solve()
<<<3.5484>>>