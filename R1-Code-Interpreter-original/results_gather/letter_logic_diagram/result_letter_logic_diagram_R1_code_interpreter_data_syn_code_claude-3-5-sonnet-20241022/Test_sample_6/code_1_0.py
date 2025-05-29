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
        # Get the letter that should be on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] not in [None, '']:
                diag_letter = grid[i][6-i]
                break
        if diag_letter is None:
            diag_letter = 'a'  # We'll use 'a' for the diagonal
        if letter != diag_letter:
            return False
    
    return True

def solve_grid(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip if cell is already filled in initial grid
    if initial_grid[row][col] not in [None, '']:
        if col == 6:
            return solve_grid(grid, initial_grid, row + 1, 0)
        return solve_grid(grid, initial_grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            next_row = row + 1 if col == 6 else row
            next_col = 0 if col == 6 else col + 1
            if solve_grid(grid, initial_grid, next_row, next_col):
                return True
            grid[row][col] = None
    
    return False

# Initialize grid with given values
initial_grid = [
    ['c','b','a','d',None,'f',None],
    ['b',None,None,'e',None,None,None],
    [None,'d',None,'f','g',None,None],
    ['d',None,None,'g',None,None,None],
    [None,None,'g',None,None,None,None],
    [None,'g','c',None,'a',None,None],
    [None,None,'b',None,'d',None,'f']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve the puzzle
if solve_grid(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")