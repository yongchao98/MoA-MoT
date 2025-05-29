def generate_output_grid(input_array):
    # Count non-9 numbers
    non_nine_count = sum(1 for x in input_array if x != 9)
    
    # Determine size
    if non_nine_count == 1:
        size = 5
    else:
        size = 10  # Default for multiple non-9 numbers
    
    # Initialize grid with 9s
    grid = [[9 for _ in range(size)] for _ in range(size)]
    
    # Fill the grid with diagonal pattern
    for i in range(size):
        for j in range(size):
            if i == 0:  # First row
                grid[i][j] = 9
            elif j == 0:  # First column after first row
                pos = (size - i) % len(input_array)
                grid[i][j] = input_array[pos]
            else:  # Rest of the grid
                grid[i][j] = grid[i-1][j-1]
    
    # Place initial sequence in first column
    for i in range(len(input_array)):
        if input_array[i] != 9:
            grid[i][0] = input_array[i]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)