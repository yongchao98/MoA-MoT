def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    
    # Create a copy of the input grid
    output_grid = [row[:] for row in input_grid]
    
    # Process each cell
    for i in range(height):
        for j in range(width):
            if input_grid[i][j] == 0:
                # If in upper half (including middle), use number from first row
                if i <= (height-1)//2:
                    output_grid[i][j] = input_grid[0][j]
                else:
                    output_grid[i][j] = input_grid[-1][j]
    
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