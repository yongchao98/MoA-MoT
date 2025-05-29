def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_constraints(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter, initial_grid):
    # Must match initial grid
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
    
    # Check minor diagonal - must match existing diagonal letters
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    # Find first empty position
    empty = find_empty(grid)
    if not empty:
        return check_initial_constraints(grid, initial_grid)
    
    row, col = empty
    
    # Get diagonal letter if we're on diagonal
    diagonal_letter = None
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter:
            # If we have a diagonal letter, we must use it
            if is_valid(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each letter
    letters = [diagonal_letter] if diagonal_letter else 'abcdefg'
    for letter in letters:
        if letter and is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with all constraints
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# First verify that initial grid is valid
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

if valid_initial and solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")