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
    
    # Fill first row - all 9s except possibly first position
    for i in range(len(input_array)):
        if input_array[-(i+1)] != 9:
            grid[0][0] = input_array[-(i+1)]
            break
    
    # Fill subsequent rows
    for i in range(1, size):
        # Place the sequence starting from the left
        for j in range(min(i+1, len(input_array))):
            grid[i][j] = input_array[j]
            
        # Fill the rest with 9s
        for j in range(len(input_array), size):
            grid[i][j] = 9
    
    # Print the result
    result = ""
    for row in grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_array = [6, 7, 4, 5, 9]
generate_output_grid(input_array)