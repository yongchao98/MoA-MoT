def extract_subgrid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Determine the size of the output grid based on examples
    # Example 1: 10x8 -> 6x4
    # Example 2: 6x6 -> 2x2
    # Example 3: 8x10 -> 4x6
    # Example 4: 6x6 -> 2x2
    # Example 5: 10x10 -> 6x6
    
    # Calculate the number of rows and columns to remove
    # This is based on the pattern observed in the examples
    remove_rows = (input_rows - 6) // 2
    remove_cols = (input_cols - 6) // 2
    
    # Extract the subgrid
    output_grid = [row[remove_cols:input_cols-remove_cols] for row in input_grid[remove_rows:input_rows-remove_rows]]
    
    return output_grid

# Test input grid
input_grid = [
    [7, 1, 9, 9, 9, 9, 1, 6],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [9, 1, 0, 9, 0, 9, 1, 9],
    [9, 1, 9, 0, 9, 0, 1, 9],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [5, 1, 9, 9, 9, 9, 1, 8]
]

# Get the output grid
output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))