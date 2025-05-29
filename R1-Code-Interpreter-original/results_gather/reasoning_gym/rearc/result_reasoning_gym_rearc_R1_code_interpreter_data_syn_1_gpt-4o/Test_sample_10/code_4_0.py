def extract_center_subgrid(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the center row and column
    center_row = input_grid[rows // 2]
    center_col = [input_grid[i][cols // 2] for i in range(rows)]
    
    # Find unique elements in the center row and column
    unique_row_elements = list(set(center_row))
    unique_col_elements = list(set(center_col))
    
    # Determine the size of the output grid
    subgrid_size = min(len(unique_row_elements), len(unique_col_elements))
    
    # Calculate the starting index for the subgrid
    start_row = (rows - subgrid_size) // 2
    start_col = (cols - subgrid_size) // 2
    
    # Extract the subgrid
    output_grid = []
    for i in range(start_row, start_row + subgrid_size):
        output_grid.append(input_grid[i][start_col:start_col + subgrid_size])
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4],
    [4, 6, 8, 6, 4],
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4]
]

# Extract the output grid
output_grid = extract_center_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))