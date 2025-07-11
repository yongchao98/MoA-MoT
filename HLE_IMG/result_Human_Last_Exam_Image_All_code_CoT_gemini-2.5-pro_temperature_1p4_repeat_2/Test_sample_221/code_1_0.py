import math

def solve_goldberg():
    """
    This function calculates the properties of the Goldberg polyhedron shown in the image.
    """
    # Step 1 & 2: Determine m and n by visual inspection of the polyhedron's structure.
    # The path between two pentagons involves taking m steps, turning 60 degrees, and taking n steps.
    # By tracing the path on the image, we find a path of 5 steps and 3 steps.
    # Given the hint m > n, we have:
    m = 5
    n = 3

    # Step 3: Determine the number of pentagonal faces, P.
    # All Goldberg polyhedra have 12 pentagonal faces.
    P = 12

    # Step 4: Calculate the number of hexagonal faces, H.
    # The formula for the number of hexagons in a G(m,n) polyhedron is:
    # H = 10 * (m^2 + m*n + n^2 - 1)
    T = m**2 + m*n + n**2
    H = 10 * (T - 1)

    # Step 5: Format and print the output as 'm,n,H,P'.
    print(f"{m},{n},{H},{P}")

solve_goldberg()