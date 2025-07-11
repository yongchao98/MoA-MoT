def solve_grid_puzzle():
    """
    Solves the grid puzzle by applying a transformation rule.
    The rule is: the new value of a cell is 1 minus the parity of the sum
    of the cell itself and its four cardinal neighbors.
    out[r][c] = 1 - ((in[r][c] + North + South + East + West) % 2)
    """
    
    # The 6th test input grid
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(grid)
    cols = len(grid[0])
    out_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Iterate through each cell of the grid
    for r in range(rows):
        for c in range(cols):
            # Start with the cell's own value
            total_sum = grid[r][c]

            # Add cardinal neighbors' values
            # North
            if r > 0:
                total_sum += grid[r-1][c]
            # South
            if r < rows - 1:
                total_sum += grid[r+1][c]
            # West
            if c > 0:
                total_sum += grid[r][c-1]
            # East
            if c < cols - 1:
                total_sum += grid[r][c+1]
            
            # Apply the identified rule
            out_grid[r][c] = 1 - (total_sum % 2)

    # Flatten the output grid into a single string
    result_string = "".join(str(cell) for row in out_grid for cell in row)
    
    # Print the final 25-digit binary string
    print(result_string)

solve_grid_puzzle()