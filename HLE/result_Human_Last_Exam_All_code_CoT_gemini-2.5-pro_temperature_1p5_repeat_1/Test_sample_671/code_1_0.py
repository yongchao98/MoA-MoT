def print_formula():
    """
    This function prints the derived formula for the polynomial sequence f_n(p).
    The formula is constructed using f-strings to explicitly show the numerical coefficients.
    """
    
    # The numerical coefficients in the formula
    c1 = 1
    c2 = 1
    c3 = 2
    c4 = 1

    # Printing the formula in a formatted string.
    # The formula is f_n(p) = (p^n - (1 - p)^n) / (2*p - 1)
    print("The derived formula for f_n(p) is:")
    print(f"f_n(p) = (p**n - ({c1} - p)**n) / ({c3}*p - {c4})")

print_formula()