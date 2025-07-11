import math

def display_hamiltonicity_threshold():
    """
    This function explains and prints the formula for the d-threshold for Hamiltonicity.
    The problem defines d = n/2 - eta, with 1/2 <= eta <= n/64.
    The threshold p is the minimum probability such that H_n U G(n,p) is a.a.s. Hamiltonian
    for any graph H_n with minimum degree at least d.
    
    Based on results from random graph theory, the threshold p is asymptotically
    proportional to n / eta^2.
    """
    
    # Exponents in the formula p is proportional to n^(n_exponent) * eta^(eta_exponent)
    n_exponent = 1
    eta_exponent = -2
    
    print("The d-threshold for Hamiltonicity, p, is described by the following relation:")
    print("p = C * (n / (eta * eta))")
    print("where:")
    print("  - C is a constant factor.")
    print("  - n is the number of vertices.")
    print("  - eta is given by d = n/2 - eta, where d is the minimum degree.")
    print("\nThe final equation can be expressed in terms of powers of n and eta.")
    
    # As per the instruction to output each number in the final equation,
    # we print the exponents of n and eta.
    print("\nThe numbers (exponents) in the final equation are:")
    print(f"  - The exponent of n is: {n_exponent}")
    print(f"  - The exponent of eta is: {eta_exponent}")

if __name__ == '__main__':
    display_hamiltonicity_threshold()