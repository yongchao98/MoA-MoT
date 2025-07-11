import math

def solve_goldberg_polyhedron():
    """
    This function calculates the properties of the Goldberg polyhedron shown in the image.
    """
    # Step 1: Determine m and n from visual inspection.
    # By tracing the path between two adjacent pentagons (centers of the spirals) on the hexagonal grid,
    # we can count the steps. A path of 4 steps in one direction, a 60-degree turn,
    # and 3 steps in the new direction leads to the next pentagon.
    # The hint m > n is also satisfied.
    m = 4
    n = 3

    # Step 2: The number of pentagonal faces (P) in any Goldberg polyhedron is always 12.
    P = 12

    # Step 3: Calculate the number of hexagonal faces (H) using the formula.
    # The formula for H is 10 * (m^2 + m*n + n^2 - 1).
    H = 10 * (m**2 + m * n + n**2 - 1)

    # Step 4: Print the result in the required format "m,n,H,P".
    print(f"{m},{n},{H},{P}")

solve_goldberg_polyhedron()