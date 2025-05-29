def is_valid(grid, row, col, letter, initial_grid):
    # Strict check for initial positions
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row uniqueness
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column uniqueness
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Get minor diagonal letter if exists
    minor_diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diag_letter = grid[i][6-i]
            break

    # If on minor diagonal, must match existing minor diagonal letter
    if row + col == 6:
        if minor_diag_letter and letter != minor_diag_letter:
            return False

    return True

def print_solution(grid):
    for row in grid:
        print(','.join(row))

def solve(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True

    row = pos // 7
    col = pos % 7

    # If position is pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)

    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''

    return False

# Initialize grids
initial_grid = [
    ['', 'g', 'f', '', 'e', '', ''],
    ['', 'f', 'a', '', 'c', '', 'b'],
    ['', '', '', 'c', 'd', '', ''],
    ['', '', 'c', 'd', 'b', '', ''],
    ['e', 'c', '', '', '', '', 'a'],
    ['', 'd', '', 'g', 'f', 'a', 'e'],
    ['', '', 'g', 'f', '', 'e', 'c']
]

grid = [['' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    print_solution(grid)
else:
    print("No solution found")