import math

def min_vector_field_zeros(chi_M, chi_partial_M):
    """
    Calculates the least number of zeros a vector field can have on a compact
    manifold M with boundary ∂M.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.
        chi_partial_M (int): The Euler characteristic of the boundary ∂M.
    """
    # A necessary topological condition is that the Euler characteristic of the boundary
    # of a compact manifold is always an even integer.
    if chi_partial_M % 2 != 0:
        print(f"Error: The Euler characteristic of the boundary (chi_partial_M={chi_partial_M}) must be an even number.")
        return

    # The formula for the minimum number of zeros is |χ(M) - (1/2)χ(∂M)|
    term2 = chi_partial_M / 2
    result = abs(chi_M - term2)
    
    # The result must be an integer.
    final_result = int(result)

    print("The formula for the minimum number of zeros is |χ(M) - (1/2) * χ(∂M)|.")
    print("Plugging in the given values:")
    # The final code needs to output each number in the final equation.
    print(f"|{chi_M} - (1/2) * {chi_partial_M}| = |{chi_M} - {term2}| = {final_result}")
    print(f"\nThe least number of zeros is: {final_result}")

if __name__ == '__main__':
    # Example 1: M = Disk (D^2). dim=2 (even).
    # χ(M) = 1. ∂M is a circle S^1, so χ(∂M) = 0.
    print("--- Example 1: M = Disk (D^2) ---")
    min_vector_field_zeros(chi_M=1, chi_partial_M=0)
    print("\n" + "="*40 + "\n")

    # Example 2: M = Solid 3-ball (D^3). dim=3 (odd).
    # χ(M) = 1. ∂M is a sphere S^2, so χ(∂M) = 2.
    print("--- Example 2: M = Solid Ball (D^3) ---")
    min_vector_field_zeros(chi_M=1, chi_partial_M=2)
    print("\n" + "="*40 + "\n")
    
    # Example 3: M = Annulus. dim=2 (even).
    # χ(M) = 0. ∂M is two disjoint circles, so χ(∂M) = 0 + 0 = 0.
    print("--- Example 3: M = Annulus ---")
    min_vector_field_zeros(chi_M=0, chi_partial_M=0)
    print("\n" + "="*40 + "\n")

    # Example 4: M = Surface of genus 2 with one disk removed. dim=2 (even).
    # χ(M) = -3. ∂M is one circle, so χ(∂M) = 0.
    print("--- Example 4: M = Genus 2 surface with a hole ---")
    min_vector_field_zeros(chi_M=-3, chi_partial_M=0)