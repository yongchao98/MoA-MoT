def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_value(grid):
    # Returns the first non-empty value on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check minor diagonal
    if row + col == 6:  # If this is a diagonal position
        diag_value = get_diagonal_value(grid)
        if diag_value and letter != diag_value:
            return False
    
    return True

def solve_grid(grid):
    # First, determine diagonal value if any exists
    diag_value = get_diagonal_value(grid)
    if not diag_value:
        # If no diagonal value exists, we need to try each possibility
        for d in 'abcdefg':
            if try_with_diagonal(grid, d):
                return True
        return False
    else:
        # If diagonal value exists, solve with that value
        return try_with_diagonal(grid, diag_value)

def try_with_diagonal(grid, diag_value):
    # First fill diagonal with the chosen value
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid_move(grid, i, 6-i, diag_value):
                return False
            grid[i][6-i] = diag_value
    
    # Then solve rest of the grid
    return solve_remaining(grid, 0, 0)

def solve_remaining(grid, row, col):
    if col >= 7:
        row += 1
        col = 0
    if row >= 7:
        return True
        
    # Skip pre-filled cells and diagonal cells
    if grid[row][col] != '' or row + col == 6:
        return solve_remaining(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid, row, col + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

print('<<<')
if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')