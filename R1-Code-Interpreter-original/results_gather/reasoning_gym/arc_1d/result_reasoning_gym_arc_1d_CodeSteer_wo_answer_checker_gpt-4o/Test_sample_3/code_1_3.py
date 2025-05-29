def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # No non-zero numbers found
    
    # Initialize the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the first non-zero number in the output grid
    output_grid[first_non_zero_index] = input_grid[first_non_zero_index]
    
    # Find the contiguous sequence of non-zero numbers after the first one
    sequence_start = first_non_zero_index + 1
    sequence = []
    
    for i in range(sequence_start, len(input_grid)):
        if input_grid[i] != 0:
            sequence.append(input_grid[i])
        else:
            break
    
    # Shift the sequence to start immediately after the first non-zero number
    if sequence:
        output_grid[first_non_zero_index + 1:first_non_zero_index + 1 + len(sequence) - 1] = sequence[:-1]
        # Replace the last number in the sequence with '0'
        output_grid[first_non_zero_index + 1 + len(sequence) - 1] = 0
    
    return output_grid

# Test input
input_grid = [1, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")