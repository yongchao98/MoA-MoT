def solve_puzzle():
    """
    Solves the grid transformation puzzle by applying a derived rule.
    The rule is based on a cell's current value and the sum of its 8 Moore neighbors.
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

    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of 8 Moore neighbors
            moore_sum = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    
                    nr, nc = r + dr, c + dc
                    
                    if 0 <= nr < rows and 0 <= nc < cols:
                        moore_sum += grid[nr][nc]
            
            c_in = grid[r][c]
            out_val = 0
            
            # Apply the discovered transformation rule
            if c_in == 0:
                if moore_sum in [2, 3]:
                    out_val = 1
            elif c_in == 1:
                if moore_sum in [2, 3, 5]:
                    out_val = 1
            
            output_grid[r][c] = out_val
            
    # Flatten the grid and format as a binary string
    result_string = "".join([str(cell) for row in output_grid for cell in row])
    
    # Print the final result
    print(result_string)

solve_puzzle()