import math

def min_zeros_on_manifold_with_boundary(chi_M, chi_dM):
    """
    Calculates the minimum number of zeros for a vector field on a compact
    manifold M with boundary ∂M.

    Args:
        chi_M (int): The Euler characteristic of the manifold M.
        chi_dM (int): The Euler characteristic of the boundary ∂M.
    """
    # The minimum number of zeros depends on the Euler characteristic of the double of M.
    # The Euler characteristic of the double D(M) is 2*chi(M) - chi(∂M).
    # A non-vanishing vector field exists iff chi(D(M)) = 0.
    test_value = 2 * chi_M - chi_dM

    print(f"Let chi(M) = {chi_M} and chi(∂M) = {chi_dM}.")
    print("The condition for a zero-free vector field is 2*chi(M) - chi(∂M) = 0.")
    print(f"The calculation is: 2 * {chi_M} - {chi_dM} = {test_value}")

    if test_value == 0:
        min_zeros = 0
        print("\nThe condition is met. The minimum number of zeros is 0.")
    else:
        min_zeros = 1
        print("\nThe condition is not met. The minimum number of zeros is 1.")

    return min_zeros

# Example: A 3-ball with a smaller ball removed from the interior.
# M is homotopy equivalent to a 2-sphere S^2, so chi(M) = 2.
# The boundary ∂M is two disjoint 2-spheres, so chi(∂M) = 2 + 2 = 4.
chi_M_example = 2
chi_dM_example = 4

min_zeros_on_manifold_with_boundary(chi_M_example, chi_dM_example)

print("-" * 20)
# Example 2: A disk D^2
# chi(M) = 1, chi(∂M) = chi(S^1) = 0
min_zeros_on_manifold_with_boundary(1, 0)