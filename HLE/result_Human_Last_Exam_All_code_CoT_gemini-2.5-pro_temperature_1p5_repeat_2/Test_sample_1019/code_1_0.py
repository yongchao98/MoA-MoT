def display_hamiltonicity_threshold_formula():
    """
    This function presents the formula for the d-threshold for Hamiltonicity
    for a graph H_n U G(n, p) where the minimum degree of H_n is d = n/2 - eta.
    """
    
    # The formula is p = ((2 * eta) / n)^(2 * eta - 1) / n.
    # The constant numbers in this final equation are 2 and -1.
    
    base_coefficient = 2
    exponent_coefficient = 2
    exponent_constant = -1

    # Representing the variables n and eta symbolically for printing.
    n_var = "n"
    eta_var = "η"

    print("The d-threshold for Hamiltonicity (p) is given by the formula:")
    print(f"p = (({base_coefficient} * {eta_var}) / {n_var})^({exponent_coefficient} * {eta_var} + ({exponent_constant})) / {n_var}")
    print("\n")
    print("The constant numbers in this final equation are:")
    print(f"In the base of the power, the coefficient of {eta_var} is: {base_coefficient}")
    print(f"In the exponent, the coefficient of {eta_var} is: {exponent_coefficient}")
    print(f"In the exponent, the constant term is: {exponent_constant}")

if __name__ == '__main__':
    display_hamiltonicity_threshold_formula()