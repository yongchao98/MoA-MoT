def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # First check: if this position was pre-filled, it must match
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:
        minor_diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                if minor_diag_letter is None:
                    minor_diag_letter = grid[i][j]
                elif grid[i][j] != minor_diag_letter:
                    return False
        if minor_diag_letter is not None and letter != minor_diag_letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

# Create a working grid
grid = [row[:] for row in initial_grid]

# Verify initial grid is valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            temp = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid(initial_grid, i, j, temp, initial_grid):
                valid_initial = False
            initial_grid[i][j] = temp
            if not valid_initial:
                break
    if not valid_initial:
        break

if valid_initial and solve_grid(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution found")