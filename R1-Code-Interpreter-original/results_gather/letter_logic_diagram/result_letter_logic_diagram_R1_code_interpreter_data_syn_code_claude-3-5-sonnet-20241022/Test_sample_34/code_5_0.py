def is_valid(grid, row, col, letter, original):
    # Must match pre-filled cell
    if original[row][col] != '' and letter != original[row][col]:
        return False
    
    # Must be 'g' on minor diagonal
    if row + col == 6 and letter != 'g':
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
    
    # If it's a minor diagonal cell, must use 'g'
    if row + col == 6:
        if is_valid(grid, row, col, 'g', original):
            grid[row][col] = 'g'
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

if solve(grid, original):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")