def print_upper_bound_relation():
    """
    This function prints the relationship between the maximum norm (k_k,inf)
    and the covolume (V) for Hilbert modular surfaces.

    Based on the interpretation of k_k,inf as the Chern number c_1^2,
    the upper bound is derived from the Hirzebruch-Zagier formula.
    """

    # Define the variable names as strings for the equation
    k_k_inf = "k_{k,∞}"
    V = "V"
    
    # The coefficient in the upper bound relation c_1^2 <= V is 1.
    coefficient = 1

    # Print the equation representing the upper bound
    print(f"The upper bound relationship is given by the equation:")
    print(f"{k_k_inf} <= {coefficient} * {V}")

# Execute the function to print the result
print_upper_bound_relation()