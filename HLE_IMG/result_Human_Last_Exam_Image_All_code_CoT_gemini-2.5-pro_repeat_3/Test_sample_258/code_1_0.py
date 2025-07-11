import math

def solve_geometry_problem():
    """
    Calculates the length of segment CP based on the length of segment AB.
    The geometric construction reveals that CP = AB * phi, where phi is the golden ratio.
    """
    # The given length of the segment AB.
    ab_length = 18468

    # The golden ratio, phi, is (sqrt(5) + 1) / 2.
    golden_ratio = (math.sqrt(5) + 1) / 2

    # Calculate the length of the segment CP.
    cp_length = ab_length * golden_ratio

    # Print the final equation with each number, as requested.
    # The equation is CP = AB * phi.
    print(f"Based on the geometric analysis, the formula is: CP = AB * ((sqrt(5) + 1) / 2)")
    print(f"Substituting the value of AB = {ab_length}:")
    print(f"CP = {ab_length} * {golden_ratio}")
    print(f"The calculated length of CP is: {cp_length:.4f}")

solve_geometry_problem()