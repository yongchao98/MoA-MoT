def solve_grid_puzzle():
    """
    Solves the grid transformation puzzle by applying a derived rule.
    The rule is as follows:
    - A '1' cell becomes '1' if it has exactly 2 neighbors.
    - A '0' cell becomes '1' if it has 2 or 3 neighbors.
    - Otherwise, the cell becomes '0'.
    The function prints the final 25-digit binary string for the test grid.
    """
    # Test Input Grid 6
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(grid)
    cols = len(grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Iterate through each cell of the grid to calculate the new value
    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of the 8 neighbors
            neighbor_sum = 0
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if i == 0 and j == 0:
                        continue  # Skip the cell itself
                    
                    nr, nc = r + i, c + j
                    
                    # Check boundary conditions
                    if 0 <= nr < rows and 0 <= nc < cols:
                        neighbor_sum += grid[nr][nc]
            
            # Apply the transformation rule
            current_value = grid[r][c]
            new_value = 0
            if current_value == 1:
                if neighbor_sum == 2:
                    new_value = 1
            else:  # current_value == 0
                if neighbor_sum == 2 or neighbor_sum == 3:
                    new_value = 1
            
            output_grid[r][c] = new_value

    # Flatten the grid into a single string and print
    result_string = "".join(str(cell) for row in output_grid for cell in row)
    print(result_string)

solve_grid_puzzle()