def is_valid(grid, row, col, letter, original):
    # Check if we're trying to change a pre-filled cell
    if original[row][col] != '' and letter != original[row][col]:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Find the letter that should be on the minor diagonal
    diag_letter = None
    for i in range(7):
        if original[i][6-i] != '':
            diag_letter = original[i][6-i]
            break
    
    # If this is a minor diagonal cell, it must match diag_letter
    if row + col == 6:
        if diag_letter and letter != diag_letter:
            return False
        # Check all other filled minor diagonal cells
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, original):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If it's a minor diagonal cell, we must use the diagonal letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if original[i][6-i] != '':
                diag_letter = original[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, original):
                grid[row][col] = diag_letter
                if solve(grid, original):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(grid, original):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
original = [
    ['f','','','c','a','','g'],
    ['e','b','c','a','d','g',''],
    ['','c','','','g','f','e'],
    ['c','a','','','','e','b'],
    ['a','d','g','','e','','c'],
    ['d','g','','e','','c',''],
    ['','','','b','','a','']
]

grid = [row[:] for row in original]

# First, find the letter that should be on the minor diagonal
diag_letter = None
for i in range(7):
    if original[i][6-i] != '':
        diag_letter = original[i][6-i]
        break

if solve(grid, original):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")