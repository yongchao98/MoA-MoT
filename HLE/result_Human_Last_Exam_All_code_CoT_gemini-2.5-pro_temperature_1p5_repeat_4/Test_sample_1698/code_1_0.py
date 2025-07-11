def count_singular_fibers(C_squared, Ks_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface.

    The formula is derived from comparing two ways of calculating the topological
    Euler characteristic of the surface S. The conditions on the family imply
    that it is a base-point-free pencil, leading to the formula:
    N = 12 * chi - Ks_squared + 4 * g - 4

    Args:
        C_squared (float): The self-intersection number of the curve class C. This is
                         related to the genus g but not directly used in the final
                         formula if g is provided.
        Ks_squared (float): The self-intersection number of the canonical divisor of S.
        chi (float): The holomorphic Euler characteristic of S, chi(O_S).
        g (float): The genus of a general curve in the family.

    Returns:
        None. Prints the calculation and the result.
    """
    
    # The formula for the number of singular fibers (N) is:
    # N = 12 * chi - K_S^2 + 4 * g - 4
    
    # Calculate N
    N = 12 * chi - Ks_squared + 4 * g - 4
    
    # Output the result
    print("The number of singular fibers (N) is calculated using the formula:")
    print("N = 12 * chi - K_S^2 + 4 * g - 4")
    print(f"Given the values:")
    print(f"  chi = {chi}")
    print(f"  K_S^2 = {Ks_squared}")
    print(f"  g = {g}")
    print("\nThe calculation is:")
    print(f"N = 12 * {chi} - {Ks_squared} + 4 * {g} - 4")
    print(f"N = {12 * chi} - {Ks_squared} + {4 * g} - 4")
    print(f"N = {N}")

if __name__ == '__main__':
    # You can change these values to test with different surfaces and families.
    # Example: A pencil of elliptic curves (g=1) on a K3 surface.
    # For a K3 surface, chi = 2 and K_S^2 = 0.
    # For an elliptic curve, g = 1.
    # The self-intersection C^2 = 2g - 2 = 0.
    print("--- Example: Elliptic fibration on a K3 surface ---")
    count_singular_fibers(C_squared=0, Ks_squared=0, chi=2, g=1)
    
    # Example 2: A fibration by genus 2 curves on a surface with
    # chi=1, K_S^2 = 8 (like P^1 x P^1, though we need to check if such a fibration exists)
    print("\n--- Example 2: Hypothetical case ---")
    count_singular_fibers(C_squared=2, Ks_squared=8, chi=1, g=2)