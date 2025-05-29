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
    
    return True

def solve(grid, original, diag_letter):
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
    
    # If it's a minor diagonal cell, we must use diag_letter
    if row + col == 6:
        if is_valid(grid, row, col, diag_letter, original):
            grid[row][col] = diag_letter
            if solve(grid, original, diag_letter):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal cells, try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(grid, original, diag_letter):
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

# Try each possible letter for the minor diagonal
for diag_letter in 'abcdefg':
    # Check if this diagonal letter is compatible with pre-filled diagonal cells
    valid_diag = True
    for i in range(7):
        if original[i][6-i] != '' and original[i][6-i] != diag_letter:
            valid_diag = False
            break
    
    if not valid_diag:
        continue
    
    # Create a new grid for this attempt
    grid = [row[:] for row in original]
    
    # Fill all empty diagonal cells with diag_letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    # Try to solve with this diagonal letter
    if solve(grid, original, diag_letter):
        for row in grid:
            print(','.join(row))
        break