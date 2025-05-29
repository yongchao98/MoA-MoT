def find_centered_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the smallest number different from the surrounding numbers
    min_value = float('inf')
    min_position = None
    
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            if input_grid[i][j] < min_value:
                min_value = input_grid[i][j]
                min_position = (i, j)
    
    # Extract the 3x3 grid centered around the min_position
    if min_position:
        i, j = min_position
        output_grid = [row[j-1:j+2] for row in input_grid[i-1:i+2]]
        return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4],
    [4, 4, 4, 0, 4],
    [4, 4, 4, 4, 4]
]

output_grid = find_centered_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))