def get_diagonal_letter(grid):
    # We know 'b' is on the diagonal at (0,6)
    return 'b'

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

def solve_with_fixed_diagonal(grid):
    # Find empty cell (non-diagonal first)
    empty = None
    # First try non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                empty = (i, j)
                break
        if empty:
            break
    
    # If no non-diagonal empty cells, we're done
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if letter != 'b' and is_valid(grid, row, col, letter):  # Don't use 'b' in non-diagonal positions
            grid[row][col] = letter
            if solve_with_fixed_diagonal(grid):
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

# First, fill all diagonal positions with 'b'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'b'

if solve_with_fixed_diagonal(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")