def print_grid(grid):
    for row in grid:
        print(','.join(str(cell) if cell else '' for cell in row))

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] not in [None, '']:
        return letter == initial_grid[row][col]
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal (top-right to bottom-left)
    if row + col == 6:
        # Minor diagonal must be 'a'
        return letter == 'a'
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] in [None, '']:
                return i, j
    return None

def solve_sudoku(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_sudoku(grid, initial_grid):
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
grid = [[cell for cell in row] for row in initial_grid]

# Solve the puzzle
if solve_sudoku(grid, initial_grid):
    for row in grid:
        print(','.join(str(cell) if cell else '' for cell in row))
else:
    print("No solution found")