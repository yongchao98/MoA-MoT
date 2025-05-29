def create_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Output grid size is 4N where N is the larger dimension
    output_size = 4 * max(rows, cols)
    
    # Initialize output grid with 2's
    output = [[2 for _ in range(output_size)] for _ in range(output_size)]
    
    # Process each position in input grid
    for i in range(rows):
        for j in range(cols):
            val = input_grid[i][j]
            # Calculate center position in output grid
            center_i = 4 * i + 2
            center_j = 4 * j + 2
            
            # Place the number
            output[center_i][center_j] = val
            
            # Check horizontal adjacency
            if j < cols - 1 and input_grid[i][j] == input_grid[i][j+1]:
                output[center_i][center_j+1] = val
                output[center_i][center_j+2] = val
                output[center_i][center_j+3] = val
                output[center_i][center_j+4] = val
            
            # Check vertical adjacency
            if i < rows - 1 and input_grid[i][j] == input_grid[i+1][j]:
                output[center_i+1][center_j] = val
                output[center_i+2][center_j] = val
                output[center_i+3][center_j] = val
                output[center_i+4][center_j] = val

    # Trim the output grid to remove unnecessary rows/columns of 2's
    # Find the bounds of non-2 values
    min_row = output_size
    max_row = 0
    min_col = output_size
    max_col = 0
    
    for i in range(output_size):
        for j in range(output_size):
            if output[i][j] != 2:
                min_row = min(min_row, i)
                max_row = max(max_row, i)
                min_col = min(min_col, j)
                max_col = max(max_col, j)
    
    # Extract the relevant portion
    trimmed = [row[max(0, min_col-1):max_col+2] for row in output[max(0, min_row-1):max_row+2]]
    
    # Format output as string
    result = ""
    for row in trimmed:
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