import math

def print_hamiltonicity_threshold_formula():
    """
    This function presents the formula for the d-threshold for Hamiltonicity.
    
    The problem considers a graph H_n with minimum degree d >= n/2 - eta, where
    1/2 <= eta <= n/64. The threshold is the probability p such that H_n U G(n,p)
    is asymptotically almost surely Hamiltonian.
    
    Based on recent literature in probabilistic combinatorics, the main term of
    the threshold across a wide part of the specified range for eta is given by:
    p ~ log(n/eta) / n.
    
    This function will print the formula in a structured way to satisfy the prompt's
    request to output each number in the equation.
    """
    
    log_n_coeff = 1
    log_eta_coeff = -1
    n_denom_coeff = 1
    
    print("The d-threshold p for Hamiltonicity is approximated by the formula:")
    # Printing the equation by constructing it from its constituent numbers (coefficients)
    print(f"p = ({log_n_coeff} * log(n) + ({log_eta_coeff}) * log(eta)) / ({n_denom_coeff} * n)")

# Execute the function to print the formula
print_hamiltonicity_threshold_formula()