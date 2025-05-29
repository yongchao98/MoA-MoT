def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_partial(grid, row, col):
    # Check row uniqueness
    row_values = [x for x in grid[row] if x != '']
    if len(row_values) != len(set(row_values)):
        return False
        
    # Check column uniqueness
    col_values = [grid[i][col] for i in range(7) if grid[i][col] != '']
    if len(col_values) != len(set(col_values)):
        return False
        
    # Check diagonal consistency
    diag_values = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']
    if len(diag_values) > 1 and len(set(diag_values)) > 1:
        return False
        
    return True

def solve_recursive(grid, row, col, original):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # If this position is pre-filled, move to next position
    if original[row][col] != '':
        grid[row][col] = original[row][col]
        if not is_valid_partial(grid, row, col):
            return False
        return solve_recursive(grid, row, col + 1, original)

    # If on diagonal, must match other diagonal values
    if row + col == 6:
        diag_values = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']
        if diag_values:
            letters = [diag_values[0]]
        else:
            letters = ['c']  # If no diagonal value yet, try 'c'
    else:
        letters = 'abcdefg'

    for letter in letters:
        # Skip if letter already in row or column
        if letter in grid[row] or letter in [grid[i][col] for i in range(7)]:
            continue
            
        grid[row][col] = letter
        if is_valid_partial(grid, row, col):
            if solve_recursive(grid, row, col + 1, original):
                return True
        grid[row][col] = ''

    return False

# Original puzzle
original = [
    ['', 'b', 'a', '', '', '', ''],
    ['b', '', '', 'd', '', '', ''],
    ['', '', 'd', 'f', '', '', 'b'],
    ['g', 'd', '', '', 'e', '', ''],
    ['d', 'f', 'c', 'e', 'b', 'a', ''],
    ['', '', '', '', '', '', 'd'],
    ['', 'e', '', '', '', '', 'f']
]

# Verify original positions are valid
def verify_original():
    grid = [[''] * 7 for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if original[i][j] != '':
                grid[i][j] = original[i][j]
                if not is_valid_partial(grid, i, j):
                    return False
    return True

# Create empty grid
grid = [[''] * 7 for _ in range(7)]

# First verify original positions are valid
if verify_original():
    if solve_recursive(grid, 0, 0, original):
        # Verify solution matches original positions
        valid = True
        for i in range(7):
            for j in range(7):
                if original[i][j] != '' and grid[i][j] != original[i][j]:
                    valid = False
                    break
        if valid:
            for row in grid:
                print(','.join(row))
        else:
            print("Solution doesn't match original positions")
    else:
        print("No solution found")
else:
    print("Original positions are invalid")