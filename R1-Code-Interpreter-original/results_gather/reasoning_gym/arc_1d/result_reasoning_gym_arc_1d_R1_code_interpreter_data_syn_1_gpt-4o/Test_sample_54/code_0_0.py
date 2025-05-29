def transform_grid(input_grid):
    # Find the start and end of the non-zero sequence
    start = next((i for i, x in enumerate(input_grid) if x != 0), None)
    end = next((i for i, x in enumerate(input_grid[::-1]) if x != 0), None)
    end = len(input_grid) - end - 1 if end is not None else None
    
    if start is None or end is None:
        return input_grid  # No non-zero sequence found
    
    # Length of the non-zero sequence
    length = end - start + 1
    
    # Shift positions
    shift_positions = length - 1
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence in the new position
    for i in range(length):
        output_grid[start + shift_positions + i] = input_grid[start + i]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)