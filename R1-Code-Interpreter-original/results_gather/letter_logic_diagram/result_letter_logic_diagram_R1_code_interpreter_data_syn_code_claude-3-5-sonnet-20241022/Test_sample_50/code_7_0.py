def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal requirement
    if row + col == 6:  # if on diagonal
        return letter == diag_letter
    
    return True

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # Skip pre-filled cells
    if grid[row][col]:
        return solve(grid, row, col + 1)

    # If we're on the diagonal, we must use the diagonal letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, row, col + 1):
                return True
        grid[row][col] = ''
        return False

    # Try each letter
    for letter in 'abcdefg':
        if letter != diagonal_letter and is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    # Make a copy of the initial grid
    grid = [row[:] for row in initial]
    
    # Check if this diagonal letter conflicts with pre-filled cells
    valid = True
    for i in range(7):
        if grid[i][6-i]:
            if grid[i][6-i] != diagonal_letter:
                valid = False
                break
    
    if not valid:
        continue
    
    # Pre-fill all diagonal cells with this letter
    for i in range(7):
        if not grid[i][6-i]:
            grid[i][6-i] = diagonal_letter
    
    if solve(grid):
        print_solution(grid)
        break
else:
    print("No solution exists")