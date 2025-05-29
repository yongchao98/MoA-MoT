def transform_grid(input_grid):
    # Find the start and end of the non-zero sequence
    start = next(i for i, x in enumerate(input_grid) if x != 0)
    end = len(input_grid) - next(i for i, x in enumerate(reversed(input_grid)) if x != 0)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Shift the sequence left by 4 positions
    for i in range(start, end):
        output_grid[i - 4] = input_grid[i]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
output_grid = transform_grid(input_grid)
print(output_grid)