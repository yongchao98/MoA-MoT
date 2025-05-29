def verify_prefilled(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def get_diagonal_letter(grid):
    # Find the first non-empty cell on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, initial_grid, row, col, letter):
    # Check pre-filled constraint
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
        diag_letter = get_diagonal_letter(grid)
        if diag_letter and diag_letter != letter:
            return False
    
    return True

def solve_grid(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return verify_prefilled(grid, initial_grid)
    
    # If pre-filled, skip to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve_grid(grid, initial_grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with pre-filled values
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Solve the puzzle
if solve_grid(grid, initial_grid):
    # Verify solution before printing
    if verify_prefilled(grid, initial_grid):
        print("<<<")
        for row in grid:
            print(",".join(row))
        print(">>>")
    else:
        print("Solution does not match pre-filled constraints")
else:
    print("No solution exists")