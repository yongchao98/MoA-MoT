def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # If this is a pre-filled cell, only allow the initial value
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
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def get_valid_letters(grid, row, col, used_in_row, used_in_col, initial_grid):
    if initial_grid[row][col]:
        return [initial_grid[row][col]]
    
    valid = []
    for letter in 'abcdefg':
        if (letter not in used_in_row and 
            letter not in used_in_col and 
            is_valid(grid, row, col, letter, initial_grid)):
            valid.append(letter)
    return valid

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Get letters already used in this row and column
    used_in_row = set(cell for cell in grid[row] if cell)
    used_in_col = set(grid[i][col] for i in range(7) if grid[i][col])
    
    # If on minor diagonal, must be 'e'
    if row + col == 6:
        valid_letters = ['e'] if 'e' not in used_in_row and 'e' not in used_in_col else []
    else:
        valid_letters = get_valid_letters(grid, row, col, used_in_row, used_in_col, initial_grid)
    
    for letter in valid_letters:
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
    for i in range(7):
        if grid[i][6-i] != 'e':
            return False
    
    return True

# Initialize grid
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
    if verify_solution(grid):
        print_grid(grid)
    else:
        print("Solution found but verification failed")
else:
    print("No solution exists")