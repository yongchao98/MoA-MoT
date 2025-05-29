def extract_central_block(input_grid):
    # Find the unique numbers in the grid
    unique_numbers = set()
    for row in input_grid:
        unique_numbers.update(row)
    
    # The central block number is the one that appears most frequently in the center
    # Assuming the grid is large enough and the central block is surrounded by a different number
    central_number = None
    for number in unique_numbers:
        if all(number in row for row in input_grid):
            central_number = number
            break
    
    # Extract the central block
    output_grid = []
    for row in input_grid:
        if central_number in row:
            start_index = row.index(central_number)
            end_index = len(row) - row[::-1].index(central_number)
            output_grid.append(row[start_index:end_index])
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

output_grid = extract_central_block(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))