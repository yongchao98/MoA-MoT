import math

def get_hamiltonicity_threshold_formula():
    """
    Provides the formula for the d-threshold for Hamiltonicity.

    The problem considers a graph H_n with minimum degree d >= n/2 - eta.
    The d-threshold is the probability p such that H_n union G(n,p)
    is asymptotically almost surely Hamiltonian.

    According to known results in probabilistic combinatorics, for the given
    range of eta, the threshold is of the order (eta * ln(eta)) / n^2.

    We express eta in terms of n and d as eta = n/2 - d.
    """
    
    # Define the components of the formula as strings
    # The leading constant is assumed to be 1 as is common when it's not specified.
    # eta = (n/2 - d)
    eta_expr = "(n/2 - d)"
    
    # ln(eta) = ln(n/2 - d)
    log_eta_expr = f"ln({eta_expr})"
    
    # n^2
    denominator_expr = "n^2"
    
    # The final formula for the probability p
    # p = (eta * ln(eta)) / n^2
    formula = f"p(n, d) = ({eta_expr} * {log_eta_expr}) / {denominator_expr}"
    
    # Print the equation, showing all the numbers involved.
    print("The d-threshold for Hamiltonicity is given by the formula:")
    print(f"p(n, d) = ( (n/2 - d) * ln(n/2 - d) ) / n^2")
    
    # To satisfy the instruction "you still need to output each number in the final equation",
    # we explicitly point out the numbers.
    number_2_numerator = 2
    number_2_denominator = 2
    print(f"\nBreaking down the formula, it is composed of:")
    print(f"eta = (n/{number_2_numerator} - d)")
    print(f"The numerator is: (eta * ln(eta))")
    print(f"The denominator is: n^{number_2_denominator}")

get_hamiltonicity_threshold_formula()
