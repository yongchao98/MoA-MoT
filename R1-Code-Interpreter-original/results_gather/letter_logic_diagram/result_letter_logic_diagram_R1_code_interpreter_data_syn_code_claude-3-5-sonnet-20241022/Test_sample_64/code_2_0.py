def is_valid(grid, row, col, letter, initial_grid):
    # First check if this position was pre-filled
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
    
    # Check minor diagonal consistency
    if row + col == 6:
        minor_diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                minor_diagonal_letter = grid[i][6-i]
                break
        if minor_diagonal_letter and letter != minor_diagonal_letter:
            return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve_grid(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If this position was pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, pos + 1)
    
    # If on minor diagonal, try to use existing diagonal letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            if is_valid(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve_grid(grid, initial_grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with the given puzzle
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve_grid(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")