def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules. This is equivalent to counting
    the number of independent sets on a grid graph.
    
    A brute-force approach is used by checking every possible assignment.
    """
    n = 4
    m = 4
    num_cells = n * m
    total_valid_assignments = 0

    # Iterate through all 2^(n*m) possible assignments
    for i in range(1 << num_cells):
        grid = []
        # Create the grid from the bits of integer i
        for r in range(n):
            row = []
            for c in range(m):
                # Map the k-th bit to grid cell (r, c)
                k = r * m + c
                if (i >> k) & 1:
                    row.append(1)
                else:
                    row.append(0)
            grid.append(row)

        # Check if the current assignment is valid
        is_valid = True
        for r in range(n):
            for c in range(m):
                if grid[r][c] == 1:
                    # Check right neighbor
                    if c + 1 < m and grid[r][c + 1] == 1:
                        is_valid = False
                        break
                    # Check bottom neighbor
                    if r + 1 < n and grid[r + 1][c] == 1:
                        is_valid = False
                        break
            if not is_valid:
                break
        
        if is_valid:
            total_valid_assignments += 1
    
    # The question requests printing numbers from an equation, but for a brute-force
    # search, the most relevant "equation" is the final sum itself.
    print(f"The total number of satisfying 0/1 assignments for a {n}x{m} grid is:")
    print(f"Total = {total_valid_assignments}")


solve_grid_assignments()