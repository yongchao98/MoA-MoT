import math

def min_zeros_on_manifold_with_boundary(chi_M, chi_partial_M):
    """
    Calculates the least number of zeros a vector field can have on a compact
    manifold M with a non-empty boundary.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.
        chi_partial_M (int): The Euler characteristic of the boundary of M.
    """
    # For any compact manifold M with boundary, chi(∂M) is always an even integer.
    # So the result of the formula is always an integer.
    if chi_partial_M % 2 != 0:
        print("Warning: The Euler characteristic of a boundary, chi(∂M), must be an even integer.")
        # We proceed anyway, but note the topological impossibility of such input.

    # The formula for the minimum number of zeros is |χ(M) - (1/2)χ(∂M)|
    result = abs(chi_M - 0.5 * chi_partial_M)
    
    # The result should be an integer.
    result = int(result)
    
    print("The formula for the minimum number of zeros is |χ(M) - (1/2)χ(∂M)|")
    print(f"For χ(M) = {chi_M} and χ(∂M) = {chi_partial_M}:")
    print(f"|{chi_M} - (1/2) * {chi_partial_M}| = {result}")
    print(f"The least number of zeros is: {result}")
    print("-" * 20)

if __name__ == '__main__':
    # Example 1: M is a 2D disk (like a circle filled in).
    # This is an even-dimensional manifold.
    # χ(M) = 1. The boundary is a circle S^1, so χ(∂M) = 0.
    # Expected result: |1 - (1/2)*0| = 1.
    print("Example 1: M = 2-Disk (D^2)")
    min_zeros_on_manifold_with_boundary(chi_M=1, chi_partial_M=0)

    # Example 2: M is a 3D ball (like a sphere filled in).
    # This is an odd-dimensional manifold.
    # χ(M) = 1. The boundary is a sphere S^2, so χ(∂M) = 2.
    # Expected result: |1 - (1/2)*2| = |1 - 1| = 0.
    print("Example 2: M = 3-Ball (B^3)")
    min_zeros_on_manifold_with_boundary(chi_M=1, chi_partial_M=2)

    # Example 3: M is a cylinder (S^1 x [0,1]).
    # This is an even-dimensional manifold.
    # χ(M) = 0. The boundary consists of two circles, so χ(∂M) = 0 + 0 = 0.
    # Expected result: |0 - (1/2)*0| = 0.
    print("Example 3: M = Cylinder (S^1 x I)")
    min_zeros_on_manifold_with_boundary(chi_M=0, chi_partial_M=0)
    
    # Example 4: M is a torus with a disk removed (punctured torus).
    # This is an even-dimensional manifold.
    # χ(M) = -1. The boundary is a circle S^1, so χ(∂M) = 0.
    # Expected result: |-1 - (1/2)*0| = 1.
    print("Example 4: M = Punctured Torus")
    min_zeros_on_manifold_with_boundary(chi_M=-1, chi_partial_M=0)