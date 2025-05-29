def transform_grid(input_grid):
    # Step 1: Extract non-zero numbers
    non_zero_numbers = [num for num in input_grid if num != 0]
    
    # Step 2: Calculate the starting position for non-zero numbers
    total_length = len(input_grid)
    non_zero_length = len(non_zero_numbers)
    start_position = (total_length - non_zero_length) // 2
    
    # Step 3: Create the output grid
    output_grid = [0] * total_length
    
    # Step 4: Place non-zero numbers in the calculated positions
    for i, num in enumerate(non_zero_numbers):
        output_grid[start_position + i] = num
    
    return output_grid

# Test input
input_grid = [6, 3, 4, 2, 6, 9, 8, 1, 7, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")