def transform_grid(input_grid):
    # Find the first sequence of non-zero numbers
    start_index = 0
    while start_index < len(input_grid) and input_grid[start_index] == 0:
        start_index += 1
    
    end_index = start_index
    while end_index < len(input_grid) and input_grid[end_index] != 0:
        end_index += 1
    
    # Extract the sequence
    non_zero_sequence = input_grid[start_index:end_index]
    sequence_length = len(non_zero_sequence)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[sequence_length:sequence_length + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 8, 9, 4, 2, 5, 4, 7, 7, 2, 9, 7, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")