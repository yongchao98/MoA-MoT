def solve_grid_transformation():
    """
    Solves the grid transformation puzzle by applying a cellular automaton rule.
    The rule is: a cell's next state is 1 if the sum of its 8 neighbors is 2 or 3,
    and 0 otherwise.
    """
    input_grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Iterate through each cell of the grid to calculate its new state
    for r in range(rows):
        for c in range(cols):
            neighbor_sum = 0
            # Iterate through the 8 neighbors (Moore neighborhood)
            for i in range(-1, 2):
                for j in range(-1, 2):
                    # Skip the cell itself
                    if i == 0 and j == 0:
                        continue
                    
                    nr, nc = r + i, c + j
                    
                    # Check for boundary conditions and sum up neighbors
                    if 0 <= nr < rows and 0 <= nc < cols:
                        neighbor_sum += input_grid[nr][nc]
            
            # Apply the transformation rule to determine the output cell's value
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Flatten the output grid into a single binary string for the final answer
    result_string = ""
    for row in output_grid:
        for cell in row:
            result_string += str(cell)
    
    print(result_string)

solve_grid_transformation()