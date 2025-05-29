def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_positions():
    return [
        (0, 1, 'f'), (0, 3, 'c'),
        (1, 1, 'd'), (1, 3, 'e'),
        (2, 5, 'b'),
        (3, 0, 'c'), (3, 1, 'e'), (3, 2, 'g'), (3, 4, 'b'), (3, 5, 'f'),
        (4, 4, 'f'), (4, 5, 'd'),
        (5, 1, 'a'), (5, 2, 'b'), (5, 5, 'c'),
        (6, 1, 'b'), (6, 2, 'f'), (6, 3, 'd'), (6, 6, 'g')
    ]

def is_valid(grid, row, col, letter, fixed_positions):
    # Check if this is a fixed position
    if (row, col) in [(r, c) for r, c, l in fixed_positions]:
        return letter == [l for r, c, l in fixed_positions if r == row and c == col][0]

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

def find_empty(grid, fixed_positions):
    # First try fixed positions
    for row, col, letter in fixed_positions:
        if grid[row][col] == '':
            return row, col

    # Then try minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j

    # Then try remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, fixed_positions):
    pos = find_empty(grid, fixed_positions)
    if not pos:
        return True

    row, col = pos

    # If this is a fixed position, only try the fixed letter
    if (row, col) in [(r, c) for r, c, l in fixed_positions]:
        letter = [l for r, c, l in fixed_positions if r == row and c == col][0]
        if is_valid(grid, row, col, letter, fixed_positions):
            grid[row][col] = letter
            if solve(grid, fixed_positions):
                return True
            grid[row][col] = ''
        return False

    # If this is on minor diagonal, try to use existing diagonal letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letter = grid[i][j]
                if is_valid(grid, row, col, letter, fixed_positions):
                    grid[row][col] = letter
                    if solve(grid, fixed_positions):
                        return True
                    grid[row][col] = ''
                return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, fixed_positions):
            grid[row][col] = letter
            if solve(grid, fixed_positions):
                return True
            grid[row][col] = ''

    return False

# Initialize empty grid
grid = [[''] * 7 for _ in range(7)]

# Get fixed positions
fixed_positions = get_initial_positions()

# First fill all fixed positions to verify initial configuration
initial_valid = True
for row, col, letter in fixed_positions:
    if not is_valid(grid, row, col, letter, fixed_positions):
        initial_valid = False
        break
    grid[row][col] = letter

if initial_valid and solve(grid, fixed_positions):
    print_grid(grid)
else:
    print("No solution exists")