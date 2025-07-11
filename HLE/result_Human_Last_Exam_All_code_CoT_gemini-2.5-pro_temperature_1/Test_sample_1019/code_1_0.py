def print_hamiltonicity_threshold_formula():
    """
    This function prints the d-threshold for Hamiltonicity for a graph H_n
    with minimum degree d = n/2 - eta, within the specified range for eta.
    The formula is presented in a piecewise manner, as the threshold's behavior
    depends on the growth rate of eta relative to n.
    """

    print("The d-threshold for Hamiltonicity, denoted by p, is for a graph H_n U G(n, p),")
    print("where H_n has a minimum degree d = n/2 - eta.")
    print("The given range for eta is 1/2 <= eta <= n/64.\n")
    print("The threshold is not a single formula but changes its form based on the value of eta.\n")

    # Case 1: eta is relatively small
    print("---------------------------------------------------------")
    print("Case 1: For the range 1/2 <= eta <= sqrt(n)")
    print("---------------------------------------------------------")
    print("The d-threshold p is of the order Theta(1/n).")
    print("\nBreaking down the equation p = Theta(1/n):")
    
    numerator_val_1 = 1
    denominator_expr_1 = "n"
    
    print(f"  - The numerator is the number {numerator_val_1}.")
    print(f"  - The denominator is the expression '{denominator_expr_1}'.")
    
    print("\n")
    
    # Case 2: eta is relatively large
    print("---------------------------------------------------------")
    print("Case 2: For the range sqrt(n) < eta <= n/64")
    print("---------------------------------------------------------")
    print("The d-threshold p is of the order Theta(eta^2 / n^2).")
    print("\nBreaking down the equation p = Theta(eta^2 / n^2):")
    
    numerator_base = "eta"
    numerator_exponent = 2
    denominator_base = "n"
    denominator_exponent = 2
    
    print(f"  - The numerator is '{numerator_base}' raised to the power of {numerator_exponent}.")
    print(f"  - The denominator is '{denominator_base}' raised to the power of {denominator_exponent}.")

print_hamiltonicity_threshold_formula()