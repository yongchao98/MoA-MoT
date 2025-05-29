def is_valid(grid, row, col, letter, original):
    # Check if we're trying to change a pre-filled cell
    if original[row][col] != ' ' and original[row][col] != letter:
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:  # Cell is on minor diagonal
        # Find any existing letter on minor diagonal
        for i in range(7):
            if grid[i][6-i] and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, original):
    empty = None
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                empty = (i, j)
                break
        if empty:
            break
            
    if not empty:
        return True
    
    row, col = empty
    # Sort letters to ensure deterministic behavior
    letters = sorted('abcdefg')
    
    # If on minor diagonal, use the existing diagonal letter if any exists
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                letters = [grid[i][6-i]]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(grid, original):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with the original puzzle
initial_grid = [
    ['', 'b', '', '', 'g', 'd', 'c'],
    ['', 'a', 'f', 'g', '', 'c', 'e'],
    ['a', 'f', '', '', 'c', 'e', 'b'],
    ['', 'g', 'd', '', 'e', '', ''],
    ['g', 'd', 'c', 'e', 'b', '', 'f'],
    ['', 'c', 'e', '', 'a', '', ''],
    ['', '', 'b', 'a', 'f', '', '']
]

# Create working grid
grid = [[cell if cell != '' else '' for cell in row] for row in initial_grid]

# Solve the puzzle
if solve(grid, initial_grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")