def solve_goldberg_polyhedron():
    """
    Calculates the properties of the Goldberg polyhedron based on visual inspection.
    """
    # Parameters m and n are determined by visually tracing the path
    # between two adjacent pentagons on the polyhedron's surface.
    # The path consists of m steps, a 60-degree turn, and n steps.
    # From the image, we count a path of 5 steps and 3 steps.
    # Given the constraint m > n.
    m = 5
    n = 3

    # The number of pentagonal faces in any Goldberg polyhedron is 12.
    P = 12

    # The number of hexagonal faces is calculated using the formula:
    # H = 10 * (m^2 + m*n + n^2 - 1)
    h_calculation_term = m**2 + m * n + n**2
    H = 10 * (h_calculation_term - 1)

    # Print the final answer in the format m,n,H,P without spaces.
    print(f"{m},{n},{H},{P}")

solve_goldberg_polyhedron()