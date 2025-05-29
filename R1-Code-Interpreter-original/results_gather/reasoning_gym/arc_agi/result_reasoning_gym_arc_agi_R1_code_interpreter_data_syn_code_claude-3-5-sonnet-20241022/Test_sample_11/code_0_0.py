def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    
    # Create a copy of the input grid
    output_grid = [row[:] for row in input_grid]
    
    # Get the numbers from first and last row
    top_number = input_grid[0][0]
    bottom_number = input_grid[-1][0]
    
    # Process each cell
    for i in range(1, height-1):  # Skip first and last rows
        for j in range(width):
            if input_grid[i][j] == 0:
                # If in upper half (including middle), use top number
                if i <= (height-1)//2:
                    output_grid[i][j] = top_number
                else:
                    output_grid[i][j] = bottom_number
    
    # Format the output
    result = ""
    for row in output_grid:
        result += " ".join(map(str, row)) + "\n"
    return result.strip()

# Test input
test_input = [
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,0,8,8,8,0,8,8,6],
    [4,8,8,8,0,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,0,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,8,8,0,8,8,0,6],
    [4,0,8,8,8,8,8,8,8,6]
]

print(process_grid(test_input))