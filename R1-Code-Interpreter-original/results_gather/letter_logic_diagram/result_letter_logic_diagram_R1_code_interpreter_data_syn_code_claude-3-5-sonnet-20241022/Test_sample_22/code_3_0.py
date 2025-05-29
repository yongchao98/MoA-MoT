def print_grid(grid):
    for row in grid:
        print(','.join(row))

def find_diagonal_letter(grid):
    # First find any pre-filled letter on the minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
        
    # Check if this is on minor diagonal and matches required letter
    if row + col == 6 and letter != diagonal_letter:
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
            
    return True

def solve(grid, diagonal_letter):
    # Find empty cell
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
    
    # If on minor diagonal, only try diagonal_letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
            
    return False

# Initial grid
initial_grid = [
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Create working copy
grid = [row[:] for row in initial_grid]

# Find the diagonal letter (must be 'c' based on pre-filled values)
diagonal_letter = find_diagonal_letter(grid)
if diagonal_letter is None:
    # If no diagonal letter is pre-filled, we need to find one that works
    for letter in 'abcdefg':
        # Try each letter on an empty diagonal position
        for i in range(7):
            if grid[i][6-i] == '':
                if is_valid(grid, i, 6-i, letter, letter):
                    diagonal_letter = letter
                    break
        if diagonal_letter:
            break

# Now solve with the diagonal letter constraint
if diagonal_letter and solve(grid, diagonal_letter):
    print_grid(grid)
else:
    print("No solution exists")