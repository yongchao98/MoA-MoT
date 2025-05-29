def create_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Output grid size calculation
    output_rows = 2 * rows + rows + 1
    output_cols = 2 * cols + cols + 1
    
    # Initialize output grid with 2's
    output = [[2 for _ in range(output_cols)] for _ in range(output_rows)]
    
    # Process each position in input grid
    for i in range(rows):
        for j in range(cols):
            val = input_grid[i][j]
            # Calculate center position in output grid
            center_i = 3 * i + 1
            center_j = 3 * j + 1
            
            # Place the center number
            output[center_i][center_j] = val
            
            # Check horizontal adjacency
            if j < cols - 1 and input_grid[i][j] == input_grid[i][j+1]:
                output[center_i][center_j+1] = val
                output[center_i][center_j+2] = val
                output[center_i][center_j+3] = val
            
            # Check vertical adjacency
            if i < rows - 1 and input_grid[i][j] == input_grid[i+1][j]:
                output[center_i+1][center_j] = val
                output[center_i+2][center_j] = val
                output[center_i+3][center_j] = val

    # Format output as string
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    return result.strip()

# Test input
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

print(create_output_grid(input_grid))