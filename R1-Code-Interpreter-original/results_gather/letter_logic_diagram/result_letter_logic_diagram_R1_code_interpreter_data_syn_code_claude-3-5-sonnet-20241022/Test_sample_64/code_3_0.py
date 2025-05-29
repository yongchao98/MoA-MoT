def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check if this violates initial constraints
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
            
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
                
    return True

def solve(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
        
    # Skip pre-filled cells
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, row, col + 1)
        
    # If on minor diagonal, use the existing diagonal letter if any
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            if is_valid_move(grid, row, col, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if solve(grid, initial_grid, row, col + 1):
                    return True
            grid[row][col] = ''
            return False
    
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
                return True
    grid[row][col] = ''
    return False

# Initial grid
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Create empty grid
grid = [[''] * 7 for _ in range(7)]

# First verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if not is_valid_move(initial_grid, i, j, initial_grid[i][j], initial_grid):
                valid = False
                break

if valid and solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")