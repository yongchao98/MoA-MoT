def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_constraints(grid, row, col, letter):
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
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid):
    # Find empty cell
    found = False
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                row, col = i, j
                found = True
                break
        if found:
            break
    
    if not found:
        return True
        
    # For cells on minor diagonal, use the existing letter if one exists
    if row + col == 6:
        minor_diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                minor_diag_letter = grid[i][j]
                break
        if minor_diag_letter:
            if verify_constraints(grid, row, col, minor_diag_letter):
                grid[row][col] = minor_diag_letter
                if solve(grid):
                    return True
                grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        if verify_constraints(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

grid = [row[:] for row in initial_grid]

# First, verify that initial positions are valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            if not verify_constraints(initial_grid, i, j, initial_grid[i][j]):
                valid = False
                break
    if not valid:
        break

if valid and solve(grid):
    print_grid(grid)
else:
    print("No solution exists")