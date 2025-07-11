import math

def solve_goldberg_polyhedron():
    """
    This script calculates the properties of the Goldberg polyhedron shown in the image.
    """
    
    # Step 1: Determine m and n from the image.
    # By visual inspection, we trace a path between two adjacent pentagons (marked by the colored rosettes).
    # The path consists of 5 steps along one axis, a 60-degree turn, and 3 steps along the new axis.
    # Given the hint m > n, we have:
    m = 5
    n = 3

    # Step 2: Determine the number of pentagonal faces, P.
    # For any Goldberg polyhedron, the number of pentagonal faces is always 12.
    P = 12

    # Step 3: Calculate the number of hexagonal faces, H.
    # The formula is H = 10 * (T - 1), where T is the triangulation number, T = m^2 + m*n + n^2.
    T = m**2 + m * n + n**2
    H = 10 * (T - 1)

    # Output the calculation process for H as requested.
    print(f"The parameters determined from the image are m = {m} and n = {n}.")
    print(f"The number of pentagonal faces, P, is always {P}.")
    print("The number of hexagonal faces, H, is calculated using the formula H = 10 * (m^2 + m*n + n^2 - 1):")
    print(f"H = 10 * ({m}^2 + {m}*{n} + {n}^2 - 1)")
    print(f"H = 10 * ({m**2} + {m*n} + {n**2} - 1)")
    print(f"H = 10 * ({T} - 1)")
    print(f"H = 10 * {T-1}")
    print(f"H = {H}")
    print("\n")

    # Step 4: Format the final answer as m,n,H,P
    print("The final answer in the format m,n,H,P is:")
    final_answer = f"{m},{n},{H},{P}"
    print(final_answer)

solve_goldberg_polyhedron()