def get_minor_diagonal_letter(grid, initial_grid):
    # First check if there's a pre-filled letter on the minor diagonal
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            return initial_grid[i][j]
    # If no pre-filled letter, check current grid
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def is_valid_placement(grid, row, col, letter, initial_grid):
    # Check initial constraints
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
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
        diagonal_letter = get_minor_diagonal_letter(grid, initial_grid)
        if diagonal_letter and diagonal_letter != letter:
            return False
    
    return True

def solve(grid, initial_grid):
    # Find empty position
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    # If this is a minor diagonal position
    if row + col == 6:
        diagonal_letter = get_minor_diagonal_letter(grid, initial_grid)
        if diagonal_letter:
            if is_valid_placement(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter, initial_grid):
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

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")