def match_pattern(subgrid, pattern):
    """Check if the subgrid matches the given pattern."""
    for i in range(len(pattern)):
        for j in range(len(pattern[0])):
            if subgrid[i][j] != pattern[i][j]:
                return False
    return True

def apply_transformation(grid, start_row, start_col, transformation):
    """Apply the transformation to the grid starting at the given position."""
    for i in range(len(transformation)):
        for j in range(len(transformation[0])):
            grid[start_row + i][start_col + j] = transformation[i][j]

def transform_grid(input_grid):
    # Define the pattern and transformation
    pattern = [
        [5, 5],
        [5, 5]
    ]
    transformation = [
        [6, 6],
        [6, 6]
    ]

    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]

    # Sliding window approach
    for i in range(len(input_grid) - len(pattern) + 1):
        for j in range(len(input_grid[0]) - len(pattern[0]) + 1):
            # Extract the subgrid
            subgrid = [row[j:j+len(pattern[0])] for row in input_grid[i:i+len(pattern)]]
            # Check if the subgrid matches the pattern
            if match_pattern(subgrid, pattern):
                # Apply the transformation
                apply_transformation(output_grid, i, j, transformation)

    return output_grid

# Define the input grid
input_grid = [
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 6, 6, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 6, 6, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2],
    [5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2, 5, 5, 2]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))