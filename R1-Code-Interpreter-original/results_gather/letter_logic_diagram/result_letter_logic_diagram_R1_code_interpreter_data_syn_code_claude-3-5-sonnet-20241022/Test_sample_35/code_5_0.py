def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'c'
    if row + col == 6:
        if letter != 'c':
            return False
    # If not on diagonal and letter is 'c', not allowed
    elif letter == 'c':
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve_grid(grid, positions):
    if not positions:
        return True
    
    pos = positions[0]
    row, col = pos
    
    if grid[row][col] != '':
        return solve_grid(grid, positions[1:])
    
    if row + col == 6:
        # Try only 'c' for diagonal positions
        if is_valid(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve_grid(grid, positions[1:]):
                return True
            grid[row][col] = ''
    else:
        # Try all letters except 'c' for non-diagonal positions
        for letter in 'abdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid, positions[1:]):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

# Create list of positions to fill, with diagonal positions first
positions = []
# Add diagonal positions first
for i in range(7):
    if grid[i][6-i] == '':
        positions.append((i, 6-i))
# Add remaining positions
for i in range(7):
    for j in range(7):
        if grid[i][j] == '' and i + j != 6:
            positions.append((i, j))

if solve_grid(grid, positions):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")