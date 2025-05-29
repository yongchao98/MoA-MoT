def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal - must be 'f'
    if row + col == 6 and letter != 'f':
        return False
    elif row + col != 6 and letter == 'f':
        return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # Skip pre-filled cells
    if grid[row][col]:
        return solve(grid, row, col + 1)

    # If we're on the diagonal, must use 'f'
    if row + col == 6:
        if is_valid(grid, row, col, 'f'):
            grid[row][col] = 'f'
            if solve(grid, row, col + 1):
                return True
        grid[row][col] = ''
        return False

    # Try each letter except 'f' for non-diagonal positions
    for letter in 'abcdeg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Pre-fill all diagonal positions with 'f'
grid = [row[:] for row in initial]
for i in range(7):
    if not grid[i][6-i]:
        grid[i][6-i] = 'f'

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")