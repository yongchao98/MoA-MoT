def print_upper_bound_relation():
    """
    This function derives and prints the upper bound for the maximum norm (k_k,∞)
    in relation to the covolume (V) for lattices associated with squarefree
    natural numbers N.

    The problem context is finding a short vector in an ideal lattice of a
    real quadratic number field (degree n=2). The relationship is a direct
    consequence of Minkowski's theorem on convex bodies.

    By applying the theorem with a hypercube as the convex body, we find a
    general bound for the infinity norm of a short vector x: ||x||_∞ <= V^(1/n).
    For the specified case, the degree n is 2.
    """

    # Define the symbols and numbers for the final equation
    max_norm_symbol = "k_{k,∞}"
    covolume_symbol = "V"
    exponent_numerator = 1
    exponent_denominator = 2

    # Print the final equation, showing each number as requested.
    print("The upper bound for the maximum norm in relation to the covolume is:")
    print(f"{max_norm_symbol} <= {covolume_symbol}^({exponent_numerator}/{exponent_denominator})")

# Execute the function to print the result
print_upper_bound_relation()