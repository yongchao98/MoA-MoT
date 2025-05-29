def is_valid(grid, initial_grid, row, col, letter):
    # First check if this conflicts with pre-filled cell
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
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            if i != row and grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve_grid(grid, initial_grid):
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
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, initial_grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with pre-filled values
initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

# Create working grid
grid = [row[:] for row in initial_grid]

if solve_grid(grid, initial_grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")