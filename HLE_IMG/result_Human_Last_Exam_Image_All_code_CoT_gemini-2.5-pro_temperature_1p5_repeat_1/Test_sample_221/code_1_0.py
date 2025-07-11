def solve_goldberg_polyhedron():
    """
    This function calculates the properties of the Goldberg polyhedron shown in the image.
    """
    
    # Step 1 & 2: Identify parameters m, n, and P based on visual analysis and definition.
    # By visually tracing the path between the centers of two adjacent pentagonal spirals,
    # we can determine the parameters m and n.
    # The path consists of 'm' steps straight, a 60-degree turn, and 'n' more steps.
    # From the image, m is approximately 5 and n is approximately 3. The hint m > n is satisfied.
    m = 5
    n = 3
    
    # All Goldberg polyhedra have exactly 12 pentagonal faces.
    P = 12
    
    # Step 3: Calculate the number of hexagonal faces, H.
    # The formula is H = 10 * (T - 1), where T = m^2 + m*n + n^2.
    T = m**2 + m * n + n**2
    H = 10 * (T - 1)
    
    # Step 4: Print the answer in the specified format "m,n,H,P".
    print(f"{m},{n},{H},{P}")

solve_goldberg_polyhedron()