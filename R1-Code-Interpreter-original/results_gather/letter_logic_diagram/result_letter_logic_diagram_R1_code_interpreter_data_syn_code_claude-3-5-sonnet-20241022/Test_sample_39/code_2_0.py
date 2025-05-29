def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial):
    # First check: must match initial if cell was pre-filled
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, initial, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If this is a pre-filled cell, skip to next position
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, pos + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')