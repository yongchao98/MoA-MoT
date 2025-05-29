def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check initial grid constraint
    if initial_grid[row][col] and initial_grid[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:
        # Find any existing letter on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i]:
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
        if not diag_letter and letter != 'e':  # If no letter found, must be 'e'
            return False

    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if not grid[i][6-i]:
            return (i, 6-i)
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return (i, j)
    return None

def solve(grid, initial_grid):
    pos = find_empty(grid)
    if not pos:
        return True

    row, col = pos
    letters = ['e'] if row + col == 6 else 'abcdefg'

    for letter in letters:
        if is_valid_move(grid, row, col, letter, initial_grid):
            # Make move
            grid[row][col] = letter
            
            # Try to solve with this move
            if solve(grid, initial_grid):
                return True
                
            # Undo move
            grid[row][col] = ''

    return False

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

# Create working grid and verify initial state
grid = [row[:] for row in initial_grid]

# First verify that initial grid is valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if grid[i][j]:
            # Temporarily remove the letter to check if it's valid
            letter = grid[i][j]
            grid[i][j] = ''
            if not is_valid_move(grid, i, j, letter, initial_grid):
                valid_initial = False
            grid[i][j] = letter
            if not valid_initial:
                break
    if not valid_initial:
        break

if not valid_initial:
    print("Initial grid is invalid")
else:
    # Try to solve
    if solve(grid, initial_grid):
        # Verify solution
        valid = True
        # Check rows and columns
        for i in range(7):
            row_set = set(grid[i])
            col_set = set(grid[j][i] for j in range(7))
            if len(row_set) != 7 or len(col_set) != 7:
                valid = False
                break
        
        # Check minor diagonal
        for i in range(7):
            if grid[i][6-i] != 'e':
                valid = False
                break
        
        if valid:
            print_grid(grid)
        else:
            print("Invalid solution found")
    else:
        print("No solution exists")