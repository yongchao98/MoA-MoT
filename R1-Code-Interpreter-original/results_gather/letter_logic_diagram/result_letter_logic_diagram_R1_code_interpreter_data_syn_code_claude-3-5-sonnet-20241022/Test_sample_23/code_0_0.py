def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this cell is on it
    if row + col == 6:  # Cell is on minor diagonal
        # Get the letter that should be on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] not in [None, letter]:
                diag_letter = grid[i][6-i]
                break
        if diag_letter is not None and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip if cell is already filled
    if grid[row][col] is not None:
        return solve_grid(grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = None
    
    return False

# Initialize grid with given values
initial_grid = [
    [None,'c','g','b','a','e','f'],
    [None,None,None,'a','e','f','d'],
    ['g','b',None,'e',None,None,'c'],
    ['b','a','e',None,'d','c',None],
    ['a','e','f','d','c',None,'b'],
    ['e','f','d',None,'g','b',None],
    ['f',None,None,'g',None,'a','e']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")