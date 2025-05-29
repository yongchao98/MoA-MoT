def transform_grid(input_grid):
    # Identify non-zero sequences
    sequences = []
    current_sequence = []
    
    for num in input_grid:
        if num != 0:
            current_sequence.append(num)
        else:
            if current_sequence:
                sequences.append(current_sequence)
                current_sequence = []
    
    if current_sequence:
        sequences.append(current_sequence)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place sequences at the leftmost positions
    position = 0
    for sequence in sequences:
        for num in sequence:
            output_grid[position] = num
            position += 1
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")