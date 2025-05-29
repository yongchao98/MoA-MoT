def transform_grid(input_grid):
    # Find the start and end of the non-zero sequence
    start = next(i for i, x in enumerate(input_grid) if x != 0)
    end = len(input_grid) - next(i for i, x in enumerate(reversed(input_grid)) if x != 0)
    
    # Shift the sequence 3 positions to the left
    output_grid = [0] * len(input_grid)
    output_grid[start-3:end-3] = input_grid[start:end]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)