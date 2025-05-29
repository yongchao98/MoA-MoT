def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        # Find first letter on minor diagonal
        first_diag = None
        for i in range(7):
            if grid[i][6-i]:
                first_diag = grid[i][6-i]
                break
        if first_diag and letter != first_diag:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if i + (6-i) == 6 and not grid[i][6-i]:
            return i, 6-i
    # Then fill rest
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    # If on minor diagonal and we have a letter there already, try that first
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                letters = grid[i][6-i] + ''.join(c for c in 'abcdefg' if c != grid[i][6-i])
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
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

# Solve
if solve(initial):
    # Print solution
    result = []
    for row in initial:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")