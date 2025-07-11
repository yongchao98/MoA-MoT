def solve_water_well_problem():
    """
    Calculates the time required for the water level in section 43 to begin to rise.
    The solution is based on a Breadth-First-Search (BFS) expansion of a "lake"
    starting from section 1.
    """
    # The 7x7 grid of depths
    grid = [
        [1, 5, 27, 22, 28, 40, 14],
        [39, 13, 17, 30, 41, 12, 2],
        [32, 35, 24, 25, 19, 47, 34],
        [16, 33, 10, 42, 7, 44, 18],
        [3, 8, 45, 37, 4, 21, 20],
        [15, 46, 38, 6, 26, 48, 49],
        [9, 23, 31, 29, 11, 36, 43]
    ]

    # The coordinates of section 1 (the pour point) are (0,0).
    # The coordinates of section 43 (the target) are (6,6).
    # The BFS distance (Manhattan distance) from (0,0) to (6,6) is 6+6=12.
    # This means section 43 joins the "lake" at the 12th expansion step.
    # The water level begins to rise in section 43 when the lake from all previous
    # steps (0 through 11) is completely filled.
    # We need to sum the depths of all sections (r, c) whose Manhattan distance
    # from (0,0) is <= 11. The distance is r + c.
    
    total_volume = 0
    
    rows = 7
    cols = 7
    
    # We print the equation as we calculate the sum.
    print("The total time is the volume of the lake before it spills into section 43's layer.")
    print("This is the sum of the depths of all sections (r,c) where r+c <= 11.")
    print("The calculation is:")
    
    equation_parts = []
    
    for r in range(rows):
        for c in range(cols):
            # Check if the Manhattan distance is within the limit
            if r + c <= 11:
                depth = grid[r][c]
                total_volume += depth
                equation_parts.append(str(depth))

    # Print the full equation
    print(" + ".join(equation_parts), f"= {total_volume}")
    
    print("\nAfter the sum is calculated:")
    print(f"The total time is {total_volume} minutes.")

# Execute the function to find the answer.
solve_water_well_problem()