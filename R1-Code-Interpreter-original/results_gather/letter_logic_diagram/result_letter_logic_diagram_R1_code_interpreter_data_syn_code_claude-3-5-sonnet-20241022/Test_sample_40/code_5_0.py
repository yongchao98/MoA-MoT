def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, initial_grid, row, col, letter):
    # Check if this position was pre-filled
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

def find_empty(grid):
    # First fill minor diagonal positions
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            return (i, 6-i, True)  # True indicates diagonal position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve(grid, initial_grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    
    # If it's a diagonal position, we must use the diagonal letter
    if is_diagonal:
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        
        if diag_letter and is_valid_placement(grid, initial_grid, row, col, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal positions, try all letters
    letters = 'abcdefg'
    for letter in letters:
        if is_valid_placement(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid
initial_grid = [
    ['', 'e', '', 'g', '', 'c', 'b'],
    ['', 'd', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Create working grid and copy initial values
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

# Find the letter that must appear on the minor diagonal
diag_letter = None
for i in range(7):
    if grid[i][6-i] != '':
        diag_letter = grid[i][6-i]
        break

if diag_letter:
    # Pre-fill the minor diagonal with the known letter
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] == '':
            grid[i][6-i] = diag_letter

print('<<<')
if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")
print('>>>')