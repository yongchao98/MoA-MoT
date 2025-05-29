def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:
        # Find the letter that should be on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty cell, check if solution is valid
    if not empty:
        return True
    
    row, col = empty
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
grid = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

# Solve the puzzle
if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")