def generate_output_grid(input_array):
    # Determine grid size
    non_nine_count = sum(1 for x in input_array if x != 9)
    if non_nine_count == 1:
        size = 5
    elif non_nine_count <= 3:
        size = 10
    else:
        size = 15
        
    # Initialize grid with 9s
    grid = [[9 for _ in range(size)] for _ in range(size)]
    
    # Place initial sequence in first row
    for j in range(len(input_array)):
        grid[0][j] = input_array[j]
    
    # Fill the rest of the grid
    for i in range(1, size):
        for j in range(size):
            # Calculate the source position from previous row
            prev_j = (j - 1) if j > 0 else (size - 1)
            grid[i][j] = grid[i-1][prev_j]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)