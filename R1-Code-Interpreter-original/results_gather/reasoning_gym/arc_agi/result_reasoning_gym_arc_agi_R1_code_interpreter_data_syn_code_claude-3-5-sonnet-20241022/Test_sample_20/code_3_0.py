def create_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Output grid size is exactly double the input dimensions
    output_rows = rows * 4 - 1
    output_cols = cols * 4 - 1
    
    # Initialize output grid with 2's
    output = [[2 for _ in range(output_cols)] for _ in range(output_rows)]
    
    # Process each position in input grid
    for i in range(rows):
        for j in range(cols):
            val = input_grid[i][j]
            # Calculate position in output grid
            out_i = i * 4
            out_j = j * 4
            
            # Place the number
            output[out_i][out_j] = val
            
            # Check horizontal adjacency
            if j < cols - 1 and input_grid[i][j] == input_grid[i][j+1]:
                output[out_i][out_j+1] = val
                output[out_i][out_j+2] = val
                output[out_i][out_j+3] = val
                output[out_i][out_j+4] = val
            
            # Check vertical adjacency
            if i < rows - 1 and input_grid[i][j] == input_grid[i+1][j]:
                output[out_i+1][out_j] = val
                output[out_i+2][out_j] = val
                output[out_i+3][out_j] = val
                output[out_i+4][out_j] = val

    # Trim excess rows and columns
    final_rows = rows * 2 + 1
    final_cols = cols * 2 + 1
    trimmed = [[output[i][j] for j in range(0, final_cols)] for i in range(0, final_rows)]
    
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