def is_valid(grid, row, col, letter):
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
    # First fill non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    return False

# Initial grid
grid = [
    ['c','','','','a','','b'],
    ['g','e','f','a','','','b'],
    ['','','','d','','b','g'],
    ['f','a','','b','c','g',''],
    ['','d','b','c','','e','f'],
    ['d','b','','','','f',''],
    ['b','','','','','','']
]

# Fill all diagonal positions with 'b' first
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'b'

if solve(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")