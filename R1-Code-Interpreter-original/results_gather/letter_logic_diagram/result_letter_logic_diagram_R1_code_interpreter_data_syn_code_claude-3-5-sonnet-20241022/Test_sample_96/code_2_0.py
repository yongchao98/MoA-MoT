def find_diagonal_letter(grid):
    # Check pre-filled cells on minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return 'e'  # If no pre-filled diagonal cells, use 'e'

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:
        return letter == diag_letter
    
    return True

def solve_grid(grid, diag_letter, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    if grid[row][col] != '':
        return solve_grid(grid, diag_letter, row, col + 1)
    
    # If on diagonal, only try the diagonal letter
    if row + col == 6:
        if is_valid(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            if solve_grid(grid, diag_letter, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        # Try all letters except diagonal letter if this position needs that letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diag_letter):
                grid[row][col] = letter
                if solve_grid(grid, diag_letter, row, col + 1):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['g','','','','','',''],
    ['','','c','','','b','g'],
    ['d','','','','','g',''],
    ['c','','','','','f',''],
    ['e','','b','g','f','','c'],
    ['a','b','','f','','',''],
    ['','','f','','c','','a']
]

# Find the letter that must be on the diagonal
diag_letter = find_diagonal_letter(initial_grid)

# Solve the puzzle
if solve_grid(initial_grid, diag_letter):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")