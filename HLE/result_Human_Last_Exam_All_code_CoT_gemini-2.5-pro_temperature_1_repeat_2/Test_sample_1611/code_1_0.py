def count_valid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a grid
    that satisfy the given adjacency implications.
    """
    n = 4
    m = 4
    num_cells = n * m
    valid_assignments_count = 0

    # Iterate through all 2^(n*m) possible assignments.
    for i in range(2**num_cells):
        # Create the grid from the bits of integer i.
        grid = []
        for r_idx in range(n):
            row = []
            for c_idx in range(m):
                cell_idx = r_idx * m + c_idx
                # Append 1 if the bit is set, otherwise 0.
                if (i >> cell_idx) & 1:
                    row.append(1)
                else:
                    row.append(0)
            grid.append(row)

        # Check if the current grid assignment is valid.
        is_valid = True
        for r in range(n):
            for c in range(m):
                if grid[r][c] == 1:
                    # Check horizontal implication: x_i,j -> NOT x_i,j+1
                    # This is violated if x_i,j = 1 and x_i,j+1 = 1
                    if c + 1 < m and grid[r][c+1] == 1:
                        is_valid = False
                        break
                    # Check vertical implication: x_i,j -> NOT x_i+1,j
                    # This is violated if x_i,j = 1 and x_i+1,j = 1
                    if r + 1 < n and grid[r+1][c] == 1:
                        is_valid = False
                        break
            if not is_valid:
                break
        
        if is_valid:
            valid_assignments_count += 1

    print(f"For a {n}x{m} grid, the number of different 0/1 assignments that satisfy all implications is:")
    print(valid_assignments_count)

if __name__ == '__main__':
    count_valid_assignments()