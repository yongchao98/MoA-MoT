def find_horizontal_symmetry(grid):
    n = len(grid)
    for row in range(n-1):  # Check between rows
        is_symmetrical = True
        for offset in range(min(row+1, n-row-1)):
            if grid[row-offset] != grid[row+1+offset]:
                is_symmetrical = False
                break
        if is_symmetrical and (row > 0 and row < n-1):  # Avoid edges
            return row
    return None

def generate_output(input_grid):
    n = len(input_grid)
    output = [row[:] for row in input_grid]
    
    # Find horizontal symmetry line
    sym_row = find_horizontal_symmetry(input_grid)
    
    if sym_row is not None:
        # Replace the row after symmetry line with 2's
        output[sym_row + 1] = [2] * n
    
    # Print the result
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7]
]

generate_output(test_input)