def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_constraints(grid):
    # Get all pre-filled positions and find diagonal letter
    prefilled = []
    diagonal_values = set()
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                prefilled.append((i, j, grid[i][j]))
                if i + j == 6:  # on minor diagonal
                    diagonal_values.add(grid[i][j])
    return prefilled, diagonal_values

def is_valid(grid, row, col, letter, diagonal_letter):
    # If on minor diagonal, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    return True

def solve(grid, pos, diagonal_letter):
    if pos >= 49:  # 7x7 = 49 cells
        return True

    row = pos // 7
    col = pos % 7

    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, pos + 1, diagonal_letter)

    letters = [diagonal_letter] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, pos + 1, diagonal_letter):
                return True
            grid[row][col] = ''

    return False

# Initial grid
initial_grid = [
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Create working copy
grid = [row[:] for row in initial_grid]

# Get constraints and verify initial state
prefilled, diagonal_values = get_initial_constraints(grid)

# Verify initial grid is valid
valid = True
for i, j, letter in prefilled:
    temp = grid[i][j]
    grid[i][j] = ''
    if not is_valid(grid, i, j, letter, 'c'):
        valid = False
    grid[i][j] = temp

# Determine diagonal letter (must be 'c' based on prefilled values)
diagonal_letter = 'c'

if valid and solve(grid, 0, diagonal_letter):
    print_grid(grid)
else:
    print("No solution exists")