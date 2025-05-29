def generate_output_grid(input_array):
    # Count non-9 numbers
    non_nine_count = sum(1 for x in input_array if x != 9)
    
    # Determine size
    if non_nine_count == 1:
        size = 5
    else:
        size = 10  # For this case
        
    # Initialize grid with 9s
    grid = [[9 for _ in range(size)] for _ in range(size)]
    
    # First row special case - all 9s
    grid[0] = [9] * size
    
    # Fill the rest of the grid
    for i in range(1, size):
        for j in range(size):
            if j < len(input_array):
                grid[i][j] = input_array[j]
            else:
                grid[i][j] = 9
        # Rotate the input array for next row
        input_array = input_array[-1:] + input_array[:-1]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)