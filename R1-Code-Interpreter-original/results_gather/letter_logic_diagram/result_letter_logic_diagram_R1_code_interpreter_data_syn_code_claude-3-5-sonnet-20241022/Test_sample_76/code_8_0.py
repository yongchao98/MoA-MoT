def validate_initial_grid(grid):
    # Check each pre-filled value for row and column conflicts
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                # Check row
                for k in range(7):
                    if k != j and grid[i][k] == grid[i][j]:
                        return False
                # Check column
                for k in range(7):
                    if k != i and grid[k][j] == grid[i][j]:
                        return False
    return True

def is_valid_placement(grid, row, col, letter):
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
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
                
    return True

def next_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then fill rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def get_diagonal_letter(grid):
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def solve(grid):
    pos = next_empty(grid)
    if not pos:
        return True
        
    row, col = pos
    
    # If on minor diagonal
    if row + col == 6:
        diag_letter = get_diagonal_letter(grid)
        if diag_letter:
            if is_valid_placement(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve(grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
            
    return False

# Initialize grid
grid = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

print('<<<')
if validate_initial_grid(grid) and solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No valid solution exists")
print('>>>')