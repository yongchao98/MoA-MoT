def display_polynomial_formula():
    """
    This function constructs and prints the derived closed-form formula for f_n(p).
    The formula is built as a string to clearly display all its parts,
    including the numerical constants 1 and 2 as requested.
    """
    # Define variable placeholders
    n = "n"
    p = "p"
    
    # Construct the numerator part of the formula: p^n - (1-p)^n
    numerator = f"({p}^{n} - (1 - {p})^{n})"
    
    # Construct the denominator part of the formula: 2p - 1
    denominator = f"(2*{p} - 1)"
    
    # Combine the parts into the final formula string
    formula = f"f_{n}({p}) = {numerator} / {denominator}"
    
    # Print the resulting formula
    print("The derived formula is:")
    print(formula)

# Execute the function to display the formula
display_polynomial_formula()