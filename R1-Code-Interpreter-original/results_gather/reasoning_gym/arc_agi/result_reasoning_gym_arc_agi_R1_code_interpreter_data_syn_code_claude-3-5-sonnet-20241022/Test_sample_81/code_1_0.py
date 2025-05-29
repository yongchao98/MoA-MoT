def transform_grid(input_grid):
    output_grid = []
    
    for row in input_grid:
        # Separate 6's and non-6's while preserving order
        non_sixes = [x for x in row if x != 6]
        num_sixes = row.count(6)
        
        # Create new row with 6's on left and other numbers on right
        new_row = [6] * num_sixes + non_sixes
        output_grid.append(new_row)
    
    return output_grid

# Test input
test_input = [
    [6, 2, 6, 2, 2],
    [8, 6, 6, 8, 6],
    [6, 6, 9, 6, 6],
    [1, 6, 6, 1, 6],
    [3, 6, 6, 6, 6]
]

result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))