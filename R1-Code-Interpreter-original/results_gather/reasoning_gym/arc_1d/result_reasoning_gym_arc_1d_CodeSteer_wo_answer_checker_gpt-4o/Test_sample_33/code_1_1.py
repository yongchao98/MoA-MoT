def transform_grid(input_grid):
    # Convert the input string to a list of integers
    grid = list(map(int, input_grid.split()))
    
    # Find the longest sequence of identical numbers at the start
    longest_sequence_value = grid[0]
    longest_sequence_length = 1
    
    for i in range(1, len(grid)):
        if grid[i] == longest_sequence_value:
            longest_sequence_length += 1
        else:
            break
    
    # Create the output grid by extending the sequence
    output_grid = grid[:]
    for i in range(longest_sequence_length, len(grid)):
        if grid[i] == 0:
            output_grid[i] = longest_sequence_value
        else:
            break
    
    # Convert the output list back to a string
    return ' '.join(map(str, output_grid))

# Test input
input_grid = "2 2 2 2 2 0 0 0 2 0 0 0 0 0 2 0 0 0 0"
output_grid = transform_grid(input_grid)
print(f"<<<{output_grid}>>>")