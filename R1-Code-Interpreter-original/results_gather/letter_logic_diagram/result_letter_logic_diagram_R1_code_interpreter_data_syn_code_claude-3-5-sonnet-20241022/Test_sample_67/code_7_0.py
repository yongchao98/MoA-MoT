def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Strict diagonal check - MUST be 'f' on diagonal
    if row + col == 6:  # if position is on diagonal
        return letter == 'f'
    elif letter == 'f':  # if trying to place 'f' elsewhere
        # Don't allow 'f' in non-diagonal positions
        return False
    
    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):  # Skip diagonal positions
                return (i, j)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If position is on diagonal, only try 'f'
    if row + col == 6:
        if is_valid(grid, row, col, 'f'):
            grid[row][col] = 'f'
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    else:
        # For non-diagonal positions, try all letters except 'f'
        for letter in 'abcdeg':  # Excluding 'f'
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', ''],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', '', '', 'c'],
    ['', '', '', '', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

# First verify that 'f' can be used on diagonal
can_use_f = True
for i in range(7):
    j = 6 - i
    if grid[i][j] != '' and grid[i][j] != 'f':
        can_use_f = False
        break

if can_use_f:
    # Pre-fill all known diagonal positions with 'f'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'f'
    
    if solve_grid(grid):
        print('<<<')
        for row in grid:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution exists")
else:
    print("No solution exists - diagonal constraint cannot be satisfied")