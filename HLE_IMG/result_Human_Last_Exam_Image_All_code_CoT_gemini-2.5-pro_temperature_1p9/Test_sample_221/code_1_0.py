import math

def solve_goldberg_polyhedron():
    """
    This function calculates the parameters and face counts of the Goldberg
    polyhedron shown in the image.
    """
    # Step 1: Determine m and n from visual inspection.
    # By tracing the path between two adjacent pentagons (centers of spirals),
    # we move 4 hex steps, turn, and move 2 more hex steps.
    # With the hint m > n, we have:
    m = 4
    n = 2

    # Step 2: Determine the number of pentagonal faces (P).
    # For all Goldberg polyhedra, P is constant.
    P = 12

    # Step 3: Calculate the number of hexagonal faces (H) using the formula.
    # Formula: H = 10 * (m^2 + m*n + n^2 - 1)
    m_squared = m**2
    m_times_n = m * n
    n_squared = n**2
    
    H = 10 * (m_squared + m_times_n + n_squared - 1)

    # Step 4: Present the results.
    print(f"From visual inspection of the paths between pentagons and the hint m > n, we find:")
    print(f"m = {m}")
    print(f"n = {n}")
    print("-" * 20)
    
    print(f"The number of pentagonal faces (P) is always 12.")
    print(f"P = {P}")
    print("-" * 20)

    print(f"The number of hexagonal faces (H) is calculated as:")
    print(f"H = 10 * (m^2 + m*n + n^2 - 1)")
    print(f"H = 10 * ({m}^2 + {m}*{n} + {n}^2 - 1)")
    print(f"H = 10 * ({m_squared} + {m_times_n} + {n_squared} - 1)")
    print(f"H = 10 * ({m_squared + m_times_n + n_squared} - 1)")
    print(f"H = 10 * {m_squared + m_times_n + n_squared - 1}")
    print(f"H = {H}")
    print("-" * 20)
    
    print(f"The final answer in the format m,n,H,P is:")
    print(f"{m},{n},{H},{P}")


solve_goldberg_polyhedron()