def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
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

    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_next_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_puzzle(grid, initial_grid):
    # Find empty position
    pos = find_next_empty(grid)
    if not pos:
        return True

    row, col = pos

    # Get the letter for minor diagonal if this position is on it
    minor_diag_letter = None
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                minor_diag_letter = grid[i][j]
                break

    # Try letters
    letters = ['abcdefg']
    if minor_diag_letter:
        letters = [minor_diag_letter]
    elif initial_grid[row][col] != '':
        letters = [initial_grid[row][col]]

    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, initial_grid):
            # Make move
            grid[row][col] = letter
            
            # Recursively try to solve rest of the puzzle
            if solve_puzzle(grid, initial_grid):
                return True
                
            # Undo move
            grid[row][col] = ''

    return False

# Initialize the puzzle
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

# Create working grid
grid = [row[:] for row in initial_grid]

# First, verify initial configuration is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            letter = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid_move(initial_grid, i, j, letter, initial_grid):
                valid = False
            initial_grid[i][j] = letter
            if not valid:
                break
    if not valid:
        break

if valid:
    if solve_puzzle(grid, initial_grid):
        print_grid(grid)
    else:
        print("No solution exists")
else:
    print("Initial configuration is invalid")