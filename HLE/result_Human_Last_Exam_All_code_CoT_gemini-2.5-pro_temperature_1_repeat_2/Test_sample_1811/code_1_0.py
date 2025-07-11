import math

def calculate_min_vector_field_zeros(chi_M, chi_dM):
    """
    Calculates the minimum number of zeros a vector field can have on a compact
    manifold M with boundary dM.

    Args:
        chi_M (int): The Euler characteristic of the manifold M, denoted χ(M).
        chi_dM (int): The Euler characteristic of the boundary of M, denoted χ(∂M).
                      This value must be an even integer.
    """
    if chi_dM % 2 != 0:
        print(f"Error: The Euler characteristic of a boundary (χ(∂M) = {chi_dM}) must be an even integer.")
        return

    # The minimum number of zeros is given by the formula |χ(M) - (1/2)χ(∂M)|.
    # Since chi_dM is even, we can use integer division //.
    half_chi_dM = chi_dM // 2
    min_zeros = abs(chi_M - half_chi_dM)

    # Output the final equation with all numbers.
    print(f"For a manifold M with χ(M) = {chi_M} and boundary ∂M with χ(∂M) = {chi_dM}:")
    print(f"The least number of zeros is |{chi_M} - (1/2) * {chi_dM}| = |{chi_M} - {half_chi_dM}| = {min_zeros}")
    print("-" * 20)

if __name__ == '__main__':
    # Example 1: A closed disk (M = D^2)
    # The Euler characteristic of a disk is χ(M) = 1.
    # Its boundary is a circle, ∂M = S^1, with Euler characteristic χ(∂M) = 0.
    chi_disk = 1
    chi_circle = 0
    calculate_min_vector_field_zeros(chi_disk, chi_circle)

    # Example 2: A closed interval (M = [0, 1])
    # The Euler characteristic of an interval is χ(M) = 1.
    # Its boundary is two points, ∂M = {0, 1}, with Euler characteristic χ(∂M) = 2.
    chi_interval = 1
    chi_two_points = 2
    calculate_min_vector_field_zeros(chi_interval, chi_two_points)

    # Example 3: A 3-dimensional ball (M = D^3)
    # The Euler characteristic of a 3-ball is χ(M) = 1.
    # Its boundary is a 2-sphere, ∂M = S^2, with Euler characteristic χ(∂M) = 2.
    chi_3_ball = 1
    chi_2_sphere = 2
    calculate_min_vector_field_zeros(chi_3_ball, chi_2_sphere)

    # Example 4: A cylinder (M = S^1 x [0,1])
    # The Euler characteristic of a cylinder is χ(M) = 0.
    # Its boundary is two disjoint circles, ∂M = S^1 sqcup S^1, with χ(∂M) = 0 + 0 = 0.
    chi_cylinder = 0
    chi_two_circles = 0
    calculate_min_vector_field_zeros(chi_cylinder, chi_two_circles)