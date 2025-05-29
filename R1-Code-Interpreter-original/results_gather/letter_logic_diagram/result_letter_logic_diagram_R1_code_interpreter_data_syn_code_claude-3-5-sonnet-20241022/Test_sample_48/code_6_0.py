def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check diagonal requirement
    if row + col == 6 and letter != 'g':
        return False
        
    return True

def find_next_cell(grid):
    # First priority: fill diagonal cells
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
            
    # Second priority: fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    cell = find_next_cell(grid)
    if not cell:
        return True
        
    row, col = cell
    
    # If it's a diagonal cell, only try 'g'
    if row + col == 6:
        if is_valid_move(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # Try each letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
                
    return False

# Initialize grid
grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', ''],
    ['', '', 'c', 'g', 'f', '', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', '', 'a'],
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

# Verify initial diagonal positions
for i in range(7):
    j = 6 - i
    if grid[i][j] != '' and grid[i][j] != 'g':
        print("Invalid initial grid: diagonal conflict")
        exit()

# Solve
if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")