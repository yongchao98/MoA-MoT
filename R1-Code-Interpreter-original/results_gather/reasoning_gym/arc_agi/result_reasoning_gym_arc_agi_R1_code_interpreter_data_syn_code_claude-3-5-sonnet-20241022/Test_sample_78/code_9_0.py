def generate_output_grid(input_array):
    # Count non-9 numbers
    non_nine_count = sum(1 for x in input_array if x != 9)
    
    # Determine size
    if non_nine_count == 1:
        size = 5
    else:
        size = 10
        
    # Initialize grid with 9s
    grid = [[9 for _ in range(size)] for _ in range(size)]
    
    # First row is all 9s
    grid[0] = [9] * size
    
    # Get sequence of non-9 numbers in reverse order
    sequence = [x for x in reversed(input_array) if x != 9]
    
    # Fill the diagonal pattern
    for i in range(1, size):
        for j in range(min(i+1, len(sequence))):
            grid[i][j] = sequence[j]
    
    # Fill the rest of the pattern
    for i in range(1, size):
        for j in range(i+1, min(i+5, size)):
            if j-i < len(sequence):
                grid[i][j] = sequence[j-i]
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)