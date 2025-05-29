def extract_subgrid(input_grid):
    # Define the size of the subgrid
    subgrid_size = 4
    # Get the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Start from the bottom-right corner of the grid
    for i in range(rows - 1, -1, -1):
        for j in range(cols - 1, -1, -1):
            if input_grid[i][j] != 5:
                # Determine the top-left corner of the subgrid
                start_row = max(0, i - subgrid_size + 1)
                start_col = max(0, j - subgrid_size + 1)
                # Extract the 4x4 subgrid
                subgrid = [row[start_col:start_col + subgrid_size] for row in input_grid[start_row:start_row + subgrid_size]]
                return subgrid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 2, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 0, 5, 5, 5, 5, 5, 1, 1, 5, 5],
    [5, 5, 0, 0, 5, 5, 5, 5, 5, 1, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

# Extract and print the subgrid
output_grid = extract_subgrid(input_grid)
print(output_grid)