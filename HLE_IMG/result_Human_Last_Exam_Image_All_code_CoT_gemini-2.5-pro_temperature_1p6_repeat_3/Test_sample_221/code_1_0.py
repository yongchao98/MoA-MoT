def solve_goldberg():
    """
    Calculates the properties of the Goldberg polyhedron shown in the image.
    
    The parameters m and n are determined by visual inspection of the path
    between pentagons on the polyhedron's surface.
    - Start at a pentagon.
    - Move m steps in a straight line along the hexagonal grid.
    - Turn 60 degrees.
    - Move n steps in a straight line to reach the next pentagon.
    
    From the image, we can trace a path between two adjacent pentagons,
    which corresponds to m=5 and n=1.
    """
    
    # Parameters determined by visual inspection of the image
    # We trace a path from one pentagon to the next. The path consists of 5 steps,
    # a 60-degree turn, and 1 more step.
    m = 5
    n = 1
    
    # The number of pentagonal faces in any Goldberg polyhedron is always 12.
    P = 12
    
    # The number of hexagonal faces is given by the formula:
    # H = 10 * (m^2 + m*n + n^2 - 1)
    H = 10 * (m**2 + m * n + n**2 - 1)
    
    # The required output format is m,n,H,P without spaces.
    print(f"{m},{n},{H},{P}")

solve_goldberg()