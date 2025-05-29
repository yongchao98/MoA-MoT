def transform_grid(input_grid):
    # Find the position of '2'
    pos_of_2 = input_grid.index(2)
    
    # Extract non-zero numbers except '2'
    non_zero_numbers = [num for num in input_grid if num != 0 and num != 2]
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place '2' in its original position
    output_grid[pos_of_2] = 2
    
    # Place other non-zero numbers to the leftmost available positions
    for i, num in enumerate(non_zero_numbers):
        output_grid[i] = num
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 2, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")