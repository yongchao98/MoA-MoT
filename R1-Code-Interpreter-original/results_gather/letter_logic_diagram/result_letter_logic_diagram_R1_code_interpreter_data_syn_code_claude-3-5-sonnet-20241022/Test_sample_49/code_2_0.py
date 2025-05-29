def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6:  # If on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != 'e':
                return False
        if letter != 'e':
            return False
    
    return True

def get_next_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        if grid[i][6-i] == '' and i + (6-i) == 6:
            return (i, 6-i)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    empty = get_next_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, must be 'e'
    if row + col == 6:
        if is_valid_placement(grid, row, col, 'e'):
            grid[row][col] = 'e'
            if solve(grid):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    # Optimize search by avoiding letters already in row/column
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    available = [l for l in letters if l not in used_in_row and l not in used_in_col]
    
    for letter in available:
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', '', '', '', 'd', '', 'e'],
    ['c', 'f', '', '', '', 'g', ''],
    ['f', 'b', '', '', '', 'a', 'c'],
    ['', 'd', '', 'e', '', 'c', ''],
    ['', 'g', '', 'a', '', '', 'b'],
    ['', '', 'a', 'c', '', '', ''],
    ['', '', '', '', 'b', 'd', '']
]

# Verify initial state is valid
valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            temp = grid[i][j]
            grid[i][j] = ''
            if not is_valid_placement(grid, i, j, temp):
                valid = False
            grid[i][j] = temp

if valid and solve(grid):
    print_grid(grid)
else:
    print("No solution exists")