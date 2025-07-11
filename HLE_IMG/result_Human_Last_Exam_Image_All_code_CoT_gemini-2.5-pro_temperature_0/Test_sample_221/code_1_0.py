def solve_goldberg_polyhedron():
    """
    Calculates the parameters and face counts for the given Goldberg polyhedron.
    """
    # Step 1: Define parameters m and n based on visual inspection of the image.
    # The path between adjacent pentagons consists of 4 steps, a turn, and 2 steps.
    # With the hint m > n, we have:
    m = 4
    n = 2
    print(f"From visual inspection, we determine m = {m} and n = {n}.")

    # Step 2: The number of pentagonal faces in a Goldberg polyhedron is always 12.
    P = 12
    print(f"The number of pentagonal faces (P) is always {P}.")

    # Step 3 & 4: Calculate the triangulation number (T) and then the number of hexagonal faces (H).
    # The formula for T is: T = m^2 + m*n + n^2
    T = m**2 + m * n + n**2
    print(f"Calculating the triangulation number T = {m}^2 + ({m} * {n}) + {n}^2 = {m**2} + {m*n} + {n**2} = {T}")

    # The formula for H is: H = 10 * (T - 1)
    H = 10 * (T - 1)
    print(f"Calculating the number of hexagonal faces H = 10 * (T - 1) = 10 * ({T} - 1) = {H}")

    # Step 5: Format the final answer as m,n,H,P.
    print("\nFinal Answer (m,n,H,P):")
    print(f"{m},{n},{H},{P}")

solve_goldberg_polyhedron()