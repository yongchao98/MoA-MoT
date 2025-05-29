def transform_grid(input_grid):
    # Initialize lists to hold the numbers
    nines = []
    fours = []
    fives = []
    sixes = []
    
    # Collect numbers from the input grid
    for row in input_grid:
        for num in row:
            if num == 9:
                nines.append(num)
            elif num == 4:
                fours.append(num)
            elif num == 5:
                fives.append(num)
            elif num == 6:
                sixes.append(num)
    
    # Create the output grid
    output_grid = []
    row_length = len(input_grid[0])
    
    # Add rows of 9s
    while nines:
        row = nines[:row_length]
        output_grid.append(row + [0] * (row_length - len(row)))
        nines = nines[row_length:]
    
    # Add rows of 4s
    while fours:
        row = fours[:row_length]
        output_grid.append(row + [0] * (row_length - len(row)))
        fours = fours[row_length:]
    
    # Add rows of 5s
    while fives:
        row = fives[:row_length]
        output_grid.append(row + [0] * (row_length - len(row)))
        fives = fives[row_length:]
    
    # Add rows of 6s
    while sixes:
        row = sixes[:row_length]
        output_grid.append(row + [0] * (row_length - len(row)))
        sixes = sixes[row_length:]
    
    # Fill the rest with zeros if needed
    while len(output_grid) < len(input_grid):
        output_grid.append([0] * row_length)
    
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