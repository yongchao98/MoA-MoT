def create_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Output grid size calculation
    output_size = max(rows, cols) * 2 - 1
    
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
            
            # Check vertical adjacency
            if i < rows - 1 and input_grid[i][j] == input_grid[i+1][j]:
                output[base_i+1][base_j] = val
                
            # Check for groups of same numbers
            if j < cols - 2 and input_grid[i][j] == input_grid[i][j+1] == input_grid[i][j+2]:
                output[base_i][base_j+2] = val
            if i < rows - 2 and input_grid[i][j] == input_grid[i+1][j] == input_grid[i+2][j]:
                output[base_i+2][base_j] = val

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