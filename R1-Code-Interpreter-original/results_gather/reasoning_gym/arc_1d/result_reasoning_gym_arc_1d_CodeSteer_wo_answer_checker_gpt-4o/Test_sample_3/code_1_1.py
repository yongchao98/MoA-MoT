def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # No non-zero numbers found
    
    # Collect the sequence of non-zero numbers after the first one
    non_zero_sequence = [x for x in input_grid[first_non_zero_index + 1:] if x != 0]
    
    # Create the output grid
    output_grid = input_grid[:first_non_zero_index + 1]
    
    # Add the non-zero sequence, replacing the last number with '0'
    if non_zero_sequence:
        output_grid.extend(non_zero_sequence[:-1])
        output_grid.append(0)
    
    # Fill the rest with zeros
    output_grid.extend([0] * (len(input_grid) - len(output_grid)))
    
    return output_grid

# Test input
input_grid = [1, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")