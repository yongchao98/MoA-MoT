def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    
    # Create a deep copy of the input grid
    output_grid = []
    for row in input_grid:
        output_grid.append(row[:])
    
    # For each position in the grid
    for i in range(height):
        for j in range(width):
            if input_grid[i][j] == 0:
                # If in upper half (including middle), use first row number
                if i <= 4:  # Rows 0-4 are upper half
                    output_grid[i][j] = input_grid[0][j]
                else:  # Rows 5-9 are lower half
                    output_grid[i][j] = input_grid[-1][j]
    
    return output_grid

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

result = process_grid(test_input)
for row in result:
    print(" ".join(map(str, row)))