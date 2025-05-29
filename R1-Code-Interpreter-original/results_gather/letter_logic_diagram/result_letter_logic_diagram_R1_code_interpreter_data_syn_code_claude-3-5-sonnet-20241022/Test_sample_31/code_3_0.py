def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_letter(grid):
    # Find the first non-empty cell in minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def is_valid_complete(grid, initial_grid):
    # Check all initial constraints
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def is_valid_partial(grid, row, col, letter, initial_grid):
    # Check initial constraint
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
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and diag_letter != letter:
            return False
    
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return is_valid_complete(grid, initial_grid)
    
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next cell
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)
    
    # If cell is on minor diagonal, use the existing diagonal letter if any
    if row + col == 6:
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter:
            if is_valid_partial(grid, row, col, diag_letter, initial_grid):
                grid[row][col] = diag_letter
                if solve(grid, initial_grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_partial(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['e', 'a', '', '', '', 'c', 'f'],
    ['', 'b', '', '', 'c', 'f', ''],
    ['b', 'g', 'd', 'c', 'f', '', 'a'],
    ['g', 'd', '', 'f', 'e', '', 'b'],
    ['d', '', '', 'e', '', '', 'g'],
    ['', 'f', 'e', '', 'b', '', 'd'],
    ['', '', 'a', '', 'g', 'd', 'c']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")