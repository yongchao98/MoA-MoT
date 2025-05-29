def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'b'
    if row + col == 6 and letter != 'b':
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

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    letters = ['b'] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values and force 'b' on diagonal
grid = [
    ['c','','a','e','','b','b'],  # Force last position to be 'b'
    ['','a','','f','','b',''],    # Force diagonal to be 'b'
    ['','','f','b','g','',''],
    ['','','','b','','',''],      # Force diagonal to be 'b'
    ['','','g','','b','',''],     # Force diagonal to be 'b'
    ['b','b','','d','a','e','f'], # Force second position to be 'b'
    ['b','','d','a','','f','b']   # Force first position to be 'b'
]

if solve(grid):
    print_grid(grid)
else:
    print("No solution found")