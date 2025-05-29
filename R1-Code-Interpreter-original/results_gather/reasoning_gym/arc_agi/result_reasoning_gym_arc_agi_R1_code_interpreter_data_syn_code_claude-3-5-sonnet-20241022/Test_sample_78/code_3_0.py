def generate_output_grid(input_array):
    # Count non-9 numbers
    non_nine_count = sum(1 for x in input_array if x != 9)
    
    # Determine size
    if non_nine_count == 1:
        size = 5
    elif non_nine_count <= 3:
        size = 10
    else:
        size = 15
    
    # Initialize grid with 9s
    grid = [[9 for _ in range(size)] for _ in range(size)]
    
    # Fill the grid
    for i in range(size):
        for j in range(size):
            # Calculate source position with shift
            source_pos = (j - i) % len(input_array)
            if j >= i:  # Above or on main diagonal
                grid[i][j] = input_array[source_pos]
            
    # Fill wrap-around positions
    for i in range(1, size):
        for j in range(i):
            source_pos = (j - i + len(input_array)) % len(input_array)
            grid[i][j] = input_array[source_pos]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)