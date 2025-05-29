def get_column(grid, col):
    return [grid[i][col] for i in range(7)]

def is_valid_move(grid, row, col, letter, must_be_g=False):
    # If on diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
    
    # If must_be_g is True, only 'g' is allowed
    if must_be_g and letter != 'g':
        return False
    
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in get_column(grid, col):
        return False
    
    return True

def find_empty(grid):
    # First try diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates diagonal position
    
    # Then other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    
    # If diagonal position, only try 'g'
    if is_diagonal:
        if is_valid_move(grid, row, col, 'g', True):
            grid[row][col] = 'g'
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid with initial values
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# First, verify and fill all diagonal positions with 'g'
for i in range(7):
    j = 6 - i
    if grid[i][j] != '' and grid[i][j] != 'g':
        print("Invalid initial grid: diagonal position must be 'g'")
        exit()
    if grid[i][j] == '':
        grid[i][j] = 'g'

# Now solve the rest of the grid
if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")