def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_move(grid, row, col, letter, initial):
    # Check initial constraints
    if initial[row][col] and initial[row][col] != letter:
        return False
        
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
            if grid[i][j] and grid[i][j] != letter:
                return False
                
    return True

def solve_grid(grid, initial, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next cell
    if initial[row][col]:
        grid[row][col] = initial[row][col]
        return solve_grid(grid, initial, pos + 1)
    
    # Get minor diagonal letter if exists
    minor_diag_letter = None
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j]:
                minor_diag_letter = grid[i][j]
                break
    
    # Try each possible letter
    letters = [minor_diag_letter] if minor_diag_letter else 'abcdefg'
    for letter in letters:
        if is_valid_move(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve_grid(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initialize the grid with the given puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

if solve_grid(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")