def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_constraints():
    return [
        ['', 'e', '', 'g', '', 'c', 'b'],
        ['', 'd', '', '', '', 'b', 'a'],
        ['', '', '', 'c', '', '', ''],
        ['', '', 'c', 'b', '', 'e', 'd'],
        ['f', '', '', '', 'e', 'd', ''],
        ['', '', '', '', '', '', ''],
        ['b', 'a', '', 'd', 'g', '', '']
    ]

def is_valid(grid, row, col, letter, initial):
    # Must match initial value if one exists
    if initial[row][col] != '' and initial[row][col] != letter:
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
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    # First try diagonal positions
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then try positions next to filled cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If on diagonal, try to find existing diagonal value
    if row + col == 6:
        diag_val = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_val = grid[i][6-i]
                break
        if diag_val:
            if is_valid(grid, row, col, diag_val, initial):
                grid[row][col] = diag_val
                if solve(grid, initial):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each letter in a strategic order
    # Start with letters that are already in the row/column
    used_letters = set()
    for i in range(7):
        if grid[row][i] != '':
            used_letters.add(grid[row][i])
        if grid[i][col] != '':
            used_letters.add(grid[i][col])
    
    available_letters = [l for l in 'abcdefg' if l not in used_letters]
    for letter in available_letters:
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid and copy initial values
initial = get_initial_constraints()
grid = [[''] * 7 for _ in range(7)]

# Copy initial values and verify they don't conflict
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            if not is_valid(grid, i, j, initial[i][j], initial):
                valid_initial = False
                break
            grid[i][j] = initial[i][j]
    if not valid_initial:
        break

# Find and fill diagonal letter if possible
if valid_initial:
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    if diag_letter:
        for i in range(7):
            if i + (6-i) == 6 and grid[i][6-i] == '':
                grid[i][6-i] = diag_letter

print('<<<')
if not valid_initial:
    print("No solution exists - initial configuration is invalid")
elif solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")
print('>>>')