def create_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Fixed output size for this case
    output_size = 16
    
    # Initialize output grid with 2's
    output = [[2 for _ in range(output_size)] for _ in range(output_size)]
    
    # Process each position in input grid
    for i in range(rows):
        for j in range(cols):
            val = input_grid[i][j]
            # Calculate base position in output grid
            base_i = i * 2
            base_j = j * 2
            
            # Place the number
            output[base_i][base_j] = val
            
            # Check horizontal adjacency
            if j < cols - 1 and input_grid[i][j] == input_grid[i][j+1]:
                output[base_i][base_j+1] = val
                output[base_i][base_j+2] = val
            
            # Check vertical adjacency
            if i < rows - 1 and input_grid[i][j] == input_grid[i+1][j]:
                output[base_i+1][base_j] = val
                output[base_i+2][base_j] = val
            
            # Handle special cases for groups of three
            if j < cols - 2 and input_grid[i][j] == input_grid[i][j+1] == input_grid[i][j+2]:
                output[base_i][base_j+3] = val
            if i < rows - 2 and input_grid[i][j] == input_grid[i+1][j] == input_grid[i+2][j]:
                output[base_i+3][base_j] = val

    # Format output as string
    result = ""
    for row in output[:8]:  # Only take first 8 rows for the answer
        result += " ".join(map(str, row[:16])) + "\n"  # Only take first 16 columns
    return result.strip()

# Test input
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

print(create_output_grid(input_grid))