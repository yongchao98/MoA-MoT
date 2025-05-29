def get_available_letters(grid, row, col):
    used = set()
    # Check row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    # Check column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    return [c for c in 'abcdefg' if c not in used]

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    # Check minor diagonal if applicable
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    return True

def find_next_empty(grid, initial_grid):
    # First fill minor diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '' and initial_grid[i][j] == '':
            return (i, j)
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial_grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial_grid):
    pos = find_next_empty(grid, initial_grid)
    if not pos:
        return True
    
    row, col = pos
    candidates = get_available_letters(grid, row, col)
    
    # If on minor diagonal, find existing diagonal letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
            if initial_grid[i][j] != '':
                diagonal_letter = initial_grid[i][j]
                break
        if diagonal_letter:
            candidates = [diagonal_letter] if diagonal_letter in candidates else []
    
    for letter in candidates:
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Copy initial grid
grid = [row[:] for row in initial_grid]

# Solve
if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")