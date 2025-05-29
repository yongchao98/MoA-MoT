def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Position is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def get_next_empty(grid):
    # First try to fill minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def get_available_letters(grid, row, col):
    used = set()
    # Letters used in row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    
    # Letters used in column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    # If on minor diagonal, must be same as other diagonal letters
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                return [grid[i][j]]
    
    return [l for l in 'abcdefg' if l not in used]

def solve(grid):
    pos = get_next_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    available_letters = get_available_letters(grid, row, col)
    
    for letter in available_letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'e', '', '', '', ''],
    ['g', '', 'e', 'f', '', '', 'c'],
    ['', '', '', 'a', '', '', ''],
    ['', '', '', '', '', '', 'e'],
    ['', 'a', 'd', 'c', '', '', ''],
    ['', 'd', 'c', 'g', '', '', 'b'],
    ['d', 'c', 'g', '', 'f', 'b', '']
]

# Create a deep copy to work with
working_grid = [row[:] for row in grid]

# First ensure all minor diagonal positions are 'e' (we know this from the given 'e' in row 4)
for i in range(7):
    j = 6 - i
    if working_grid[i][j] == '':
        working_grid[i][j] = 'e'
    elif working_grid[i][j] != 'e':
        print("No solution exists - diagonal conflict")
        exit()

if solve(working_grid):
    print('<<<')
    for row in working_grid:
        print(','.join(x for x in row))
    print('>>>')
else:
    print("No solution exists")