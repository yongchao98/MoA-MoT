def print_hamiltonicity_threshold():
    """
    This function prints the formula for the d-threshold for Hamiltonicity.
    """
    
    # The formula for the d-threshold p is derived from recent results in random graph theory.
    # For a graph H with minimum degree d, the threshold is p = Theta((d/n)^2 / n).
    # Substituting d = n/2 - eta gives the expression below.
    # The formula is p = (n - 2*eta)^2 / (4*n^3).
    
    # The numbers in the final equation are:
    c1 = 2
    c2 = 2
    c3 = 4
    c4 = 3
    
    # We construct and print the formula string.
    formula = f"p = (n - {c1}*eta)^{c2} / ({c3}*n^{c4})"
    
    print("The d-threshold for Hamiltonicity is given by the formula:")
    print(formula)

if __name__ == "__main__":
    print_hamiltonicity_threshold()