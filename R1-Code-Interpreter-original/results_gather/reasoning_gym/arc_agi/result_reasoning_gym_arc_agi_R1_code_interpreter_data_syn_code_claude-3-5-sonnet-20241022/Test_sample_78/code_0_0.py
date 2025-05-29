def generate_output_grid(input_array):
    # Count non-9 numbers to help determine grid size
    non_nine_count = sum(1 for x in input_array if x != 9)
    
    # Determine grid size (minimum 5, or larger if needed for pattern)
    grid_size = max(5, len(input_array) + non_nine_count)
    
    # Initialize grid with 9s
    grid = [[9 for _ in range(grid_size)] for _ in range(grid_size)]
    
    # Place numbers in diagonal pattern
    for i in range(grid_size):
        for j in range(len(input_array)):
            pos = (i + j) % grid_size
            if i < grid_size and pos < grid_size:
                if i == 0:  # First row gets input directly
                    grid[i][j] = input_array[j]
                else:  # Subsequent rows shift pattern
                    grid[i][pos] = input_array[j]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)