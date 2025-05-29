def find_symmetry_line(grid):
    n = len(grid)
    
    # Check for vertical symmetry (column)
    for col in range(n):
        is_symmetrical = True
        for i in range(n):
            for j in range(n):
                if j < col and grid[i][j] != grid[i][2*col-j]:
                    is_symmetrical = False
                    break
            if not is_symmetrical:
                break
        if is_symmetrical:
            return ('vertical', col)
    
    # Check for horizontal symmetry (row)
    for row in range(n):
        is_symmetrical = True
        for i in range(n):
            for j in range(n):
                if i < row and grid[i][j] != grid[2*row-i][j]:
                    is_symmetrical = False
                    break
            if not is_symmetrical:
                break
        if is_symmetrical:
            return ('horizontal', row)
    
    return None

def generate_output(input_grid):
    n = len(input_grid)
    output = [row[:] for row in input_grid]
    
    symmetry = find_symmetry_line(input_grid)
    if symmetry:
        direction, pos = symmetry
        if direction == 'vertical':
            for i in range(n):
                output[i][pos] = 2
        else:  # horizontal
            output[pos] = [2] * n
    
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