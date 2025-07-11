def print_polynomial_formula():
    """
    This function prints the derived formula for the polynomial sequence f_n(p).
    The formula is constructed using variables for its numerical components.
    """
    
    # Define the numbers present in the formula
    num_one = 1
    num_two = 2
    
    # The derived formula is f_n(p) = (p^n - (1 - p)^n) / (2*p - 1)
    # We construct the string representation of this formula.
    numerator = f"p^n - ({num_one} - p)^n"
    denominator = f"{num_two}*p - {num_one}"
    
    formula = f"({numerator}) / ({denominator})"
    
    print("The simple formula for f_n(p) is:")
    print(formula)

# Execute the function to print the formula
print_polynomial_formula()