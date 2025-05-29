def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    
    # Create a deep copy of the input grid
    output_grid = [[x for x in row] for row in input_grid]
    
    # Process each position
    for i in range(height):
        for j in range(width):
            if input_grid[i][j] == 0:
                # If in upper half (including middle), use first row number
                if i <= 4:  # Since grid is 10x10, rows 0-4 are upper half
                    output_grid[i][j] = input_grid[0][j]
                else:  # rows 5-9 are lower half
                    output_grid[i][j] = input_grid[9][j]
    
    # Convert to string format
    return "\n".join(" ".join(str(x) for x in row) for row in output_grid)

# Test case
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