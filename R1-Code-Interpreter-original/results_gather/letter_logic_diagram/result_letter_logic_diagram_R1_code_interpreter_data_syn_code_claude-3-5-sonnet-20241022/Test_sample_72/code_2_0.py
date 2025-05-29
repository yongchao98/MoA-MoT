def is_valid(grid, row, col, letter, initial_grid):
    # Check if we're trying to change a pre-filled cell
    if initial_grid[row][col]:
        return letter == initial_grid[row][col]
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:  # Cell is on minor diagonal
        # Get the required letter for minor diagonal (should be 'e' based on puzzle)
        diag_letter = 'e'  # We know it should be 'e' from the initial grid
        if letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

def verify_solution(grid):
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == 'e' for x in diag):
        return False
    
    return True

# Initialize grid with the pre-filled values
initial_grid = [
    ['e','c','b','','f','d','g'],
    ['c','','','f','','','e'],
    ['b','a','','','g','e','c'],
    ['','f','','','e','c',''],
    ['','d','','e','','',''],
    ['','','e','c','b','','f'],
    ['','','c','','','','']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# Solve puzzle
if solve(grid, initial_grid):
    # Verify solution before printing
    if verify_solution(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("Solution found but verification failed")
else:
    print("No solution exists")