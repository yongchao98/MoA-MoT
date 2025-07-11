def print_norm_covolume_relation():
    """
    Prints the formula for the upper bound on the maximum norm (k_k,inf)
    in relation to the covolume (V) for lattices from squarefree integers.

    This relationship is derived from Minkowski's theorem applied to
    2-dimensional lattices from real quadratic number fields.
    """
    
    # Define the components of the formula
    norm_variable = "k_{k,inf}"
    covolume_variable = "V"
    exponent_numerator = 1
    exponent_denominator = 2

    # Print the explanation and the formula
    print("Based on Minkowski's theorem for a 2-dimensional lattice, the upper bound is derived.")
    print("The formula relating the maximum norm and the covolume is:")
    
    # Output the final equation with each number explicitly shown
    print(f"{norm_variable} <= {covolume_variable}^({exponent_numerator}/{exponent_denominator})")

# Execute the function to display the result
print_norm_covolume_relation()