def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = next(i for i, x in enumerate(input_grid) if x != 0)
    last_non_zero = len(input_grid) - 1 - next(i for i, x in enumerate(reversed(input_grid)) if x != 0)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[first_non_zero] = input_grid[first_non_zero]
    output_grid[last_non_zero] = input_grid[last_non_zero]
    
    return output_grid

# Test input
input_grid = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)