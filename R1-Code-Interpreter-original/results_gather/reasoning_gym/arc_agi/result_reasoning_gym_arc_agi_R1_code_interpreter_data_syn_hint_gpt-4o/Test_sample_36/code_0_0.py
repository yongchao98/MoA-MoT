def transform_grid(input_grid):
    # Flatten the grid and count occurrences of each number
    flat_grid = [num for row in input_grid for num in row]
    count_9 = flat_grid.count(9)
    count_4 = flat_grid.count(4)
    count_5 = flat_grid.count(5)
    count_6 = flat_grid.count(6)
    
    # Create the output grid
    output_grid = []
    
    # Add rows of 9s
    while count_9 > 0:
        row = [9] * min(count_9, len(input_grid[0]))
        output_grid.append(row + [0] * (len(input_grid[0]) - len(row)))
        count_9 -= len(row)
    
    # Add rows of 4s
    while count_4 > 0:
        row = [4] * min(count_4, len(input_grid[0]))
        output_grid.append(row + [0] * (len(input_grid[0]) - len(row)))
        count_4 -= len(row)
    
    # Add rows of 5s
    while count_5 > 0:
        row = [5] * min(count_5, len(input_grid[0]))
        output_grid.append(row + [0] * (len(input_grid[0]) - len(row)))
        count_5 -= len(row)
    
    # Add rows of 6s
    while count_6 > 0:
        row = [6] * min(count_6, len(input_grid[0]))
        output_grid.append(row + [0] * (len(input_grid[0]) - len(row)))
        count_6 -= len(row)
    
    # Fill the rest with zeros if needed
    while len(output_grid) < len(input_grid):
        output_grid.append([0] * len(input_grid[0]))
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 9, 0, 0, 4, 0, 4, 0, 0, 5, 0, 0, 6, 6, 6],
    [0, 0, 9, 0, 4, 0, 4, 0, 5, 5, 5, 0, 0, 6, 0],
    [9, 9, 0, 0, 4, 0, 4, 0, 0, 5, 0, 0, 6, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 9, 9, 0, 4, 4, 0, 0, 5, 0, 5, 0, 0, 6, 0],
    [0, 9, 0, 0, 0, 0, 4, 0, 5, 0, 5, 0, 6, 6, 6],
    [9, 9, 0, 0, 4, 4, 0, 0, 5, 0, 5, 0, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))