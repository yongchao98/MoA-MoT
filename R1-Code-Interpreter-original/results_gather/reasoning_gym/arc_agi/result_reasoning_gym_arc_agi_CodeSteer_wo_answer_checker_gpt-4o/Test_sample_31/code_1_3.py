def extract_output_grid(input_grid):
    # The output grid is always 6 rows high and 8 columns wide
    output_height = 6
    output_width = 8
    
    # Determine the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Initialize variables to find the region of interest
    # We will look for the largest cluster of non-border values
    max_cluster_size = 0
    best_start_row = 0
    best_start_col = 0
    
    # Iterate over possible starting points for the output grid
    for start_row in range(rows - output_height + 1):
        for start_col in range(cols - output_width + 1):
            # Count the number of non-border values in the current subgrid
            cluster_size = 0
            for r in range(start_row, start_row + output_height):
                for c in range(start_col, start_col + output_width):
                    if input_grid[r][c] != input_grid[0][0]:  # Assuming border is uniform
                        cluster_size += 1
            
            # Update the best region if this one has a larger cluster
            if cluster_size > max_cluster_size:
                max_cluster_size = cluster_size
                best_start_row = start_row
                best_start_col = start_col
    
    # Extract the best subgrid
    output_grid = []
    for r in range(best_start_row, best_start_row + output_height):
        output_grid.append(input_grid[r][best_start_col:best_start_col + output_width])
    
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 9, 9, 9, 9, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 9, 9, 9, 9, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 3, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
]

# Get the output grid
output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))