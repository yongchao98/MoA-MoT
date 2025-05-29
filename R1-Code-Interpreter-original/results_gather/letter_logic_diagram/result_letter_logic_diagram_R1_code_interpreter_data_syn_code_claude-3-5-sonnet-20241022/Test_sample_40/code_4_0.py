def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, initial_grid, row, col, letter):
    # If this is an initial position, only allow the initial value
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

def verify_solution(grid, initial_grid):
    # Verify initial values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Verify rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Verify minor diagonal
    diag_set = set(grid[i][6-i] for i in range(7))
    if len(diag_set) != 1:
        return False
    
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return verify_solution(grid, initial_grid)
    
    row = pos // 7
    col = pos % 7
    
    # If this position is pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)
    
    for letter in 'abcdefg':
        if is_valid_placement(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
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

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# First, identify the minor diagonal letter
diag_letter = None
for i in range(7):
    if initial_grid[i][6-i] != '':
        diag_letter = initial_grid[i][6-i]
        break

# Pre-fill the minor diagonal if we know the letter
if diag_letter:
    for i in range(7):
        if i + (6-i) == 6:
            grid[i][6-i] = diag_letter

print('<<<')
if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")
print('>>>')