def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, row, col, letter, initial):
    # Check if violates initial configuration
    if initial[row][col] and initial[row][col] != letter:
        return False
    
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal consistency
    if row + col == 6:
        diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
        if diag and diag[0] != letter:
            return False
    
    return True

def solve(grid, initial, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # If cell is pre-filled, move to next position
    if initial[row][col]:
        grid[row][col] = initial[row][col]
        return solve(grid, initial, row, col + 1)
    
    # Get diagonal letter if we're on diagonal
    diag_letter = None
    if row + col == 6:
        diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
        if diag:
            diag_letter = diag[0]
    
    # Try each possible letter
    letters = [diag_letter] if diag_letter else 'abcdefg'
    for letter in letters:
        if is_valid_placement(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, row, col + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Create empty working grid
grid = [[''] * 7 for _ in range(7)]

# Solve and print
if solve(grid, initial):
    for row in grid:
        print(','.join(row))