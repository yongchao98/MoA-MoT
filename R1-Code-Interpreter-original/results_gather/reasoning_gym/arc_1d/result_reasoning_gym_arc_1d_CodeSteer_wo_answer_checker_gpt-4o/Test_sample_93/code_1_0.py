from collections import deque

def transform_grid_optimized(input_grid):
    # Use deque for efficient appending
    non_zero_sequence = deque()
    
    # Collect non-zero numbers
    for num in input_grid:
        if num != 0:
            non_zero_sequence.append(num)
    
    # Count the number of zeros
    zero_count = len(input_grid) - len(non_zero_sequence)
    
    # Create the output grid
    output_grid = list(non_zero_sequence) + [0] * zero_count
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 7, 0]

# Get the output grid
output_grid = transform_grid_optimized(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")