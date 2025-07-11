def find_probability_pm():
    """
    This function explains and prints the derived formula for the probability P_m.
    'm' is treated as a symbolic positive integer.
    """

    # The derived formula for P_m is 1 / (4*m + 1).
    # As per the instruction, we need to output each number in the final equation.
    # We define the numerical constants for clarity.
    
    numerator = 1
    coefficient_of_m = 4
    constant_term = 1
    
    # We print the final formula for P_m. Since m is a variable,
    # we display the formula as a formatted string.
    print("Based on the analysis, the probability P_m for a positive integer m is:")
    print(f"P_m = {numerator} / ({coefficient_of_m}*m + {constant_term})")

# Execute the function to get the result.
find_probability_pm()