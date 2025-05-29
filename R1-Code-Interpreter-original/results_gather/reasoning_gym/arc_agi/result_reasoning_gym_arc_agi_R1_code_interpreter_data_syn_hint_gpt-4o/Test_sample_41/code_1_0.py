def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Apply the transformation rule to rows 6 to 9
    for i in range(6, 10):
        # Rearrange the middle section of the row
        # Swap or replace numbers based on a specific pattern
        # For example, swap the first and last elements of the middle section
        middle_section = input_grid[i][6:10]
        middle_section = middle_section[::-1]  # Reverse the middle section
        output_grid[i][6:10] = middle_section
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 7, 5, 2, 7, 9, 8, 8, 9, 7, 2, 5, 7, 4, 4],
    [4, 5, 5, 7, 7, 7, 8, 8, 8, 8, 7, 7, 7, 5, 5, 4],
    [7, 5, 8, 2, 9, 8, 8, 3, 3, 3, 8, 9, 2, 8, 5, 7],
    [5, 7, 2, 5, 8, 8, 5, 3, 3, 3, 8, 8, 5, 2, 7, 5],
    [2, 7, 3, 3, 3, 3, 8, 3, 3, 3, 4, 7, 8, 9, 7, 2],
    [7, 7, 3, 3, 3, 3, 8, 0, 0, 8, 4, 4, 8, 8, 7, 7],
    [9, 8, 3, 3, 3, 3, 0, 0, 0, 0, 8, 8, 5, 8, 8, 9],
    [8, 8, 3, 3, 3, 3, 0, 4, 4, 0, 0, 8, 5, 5, 8, 8],
    [8, 8, 5, 5, 8, 0, 0, 4, 4, 0, 0, 8, 5, 5, 8, 8],
    [9, 8, 8, 5, 8, 8, 0, 0, 0, 0, 8, 8, 5, 8, 8, 9],
    [7, 7, 8, 8, 4, 4, 8, 0, 0, 8, 4, 4, 8, 8, 7, 7],
    [2, 7, 9, 8, 7, 4, 8, 8, 8, 8, 4, 7, 8, 9, 7, 2],
    [5, 7, 2, 5, 8, 8, 5, 5, 5, 5, 8, 8, 5, 2, 7, 5],
    [7, 5, 8, 2, 9, 8, 8, 5, 5, 8, 8, 9, 2, 8, 5, 7],
    [4, 5, 5, 7, 7, 7, 8, 8, 8, 8, 7, 7, 7, 5, 5, 4],
    [4, 4, 7, 5, 2, 7, 9, 8, 8, 9, 7, 2, 5, 7, 4, 4]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))