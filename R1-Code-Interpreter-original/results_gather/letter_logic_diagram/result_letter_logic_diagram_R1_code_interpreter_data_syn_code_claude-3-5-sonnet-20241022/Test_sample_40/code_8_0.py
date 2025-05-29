def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this is an initial position
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
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
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def get_next_empty(grid, current_pos):
    for pos in range(current_pos, 49):
        row = pos // 7
        col = pos % 7
        if grid[row][col] == '':
            return pos
    return None

def solve(grid, initial_grid, current_pos=0):
    next_pos = get_next_empty(grid, current_pos)
    if next_pos is None:
        return True
    
    row = next_pos // 7
    col = next_pos % 7
    
    # If this is a diagonal position, find the required letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid, next_pos + 1):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, next_pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grids
initial_grid = [
    ['', 'e', '', 'g', '', 'c', 'b'],
    ['', 'd', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# Find diagonal letter
diag_letter = None
for i in range(7):
    if initial_grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

# Pre-fill diagonal if we have the letter
if diag_letter:
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            grid[i][6-i] = diag_letter

# Verify initial values are valid
initial_valid = True
for i in range(7):
    row_letters = [x for x in grid[i] if x != '']
    if len(row_letters) != len(set(row_letters)):
        initial_valid = False
        break
    col_letters = [grid[j][i] for j in range(7) if grid[j][i] != '']
    if len(col_letters) != len(set(col_letters)):
        initial_valid = False
        break

print('<<<')
if not initial_valid:
    print("No solution exists - initial configuration is invalid")
elif solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")
print('>>>')