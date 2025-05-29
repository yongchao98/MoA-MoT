def print_grid(grid):
    for row in grid:
        print(','.join(row))

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
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def get_diagonal_letter(grid):
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def solve(grid):
    # Find empty position
    empty = find_empty(grid)
    if not empty:
        return True

    row, col = empty
    diagonal_letter = get_diagonal_letter(grid)

    # If on minor diagonal and we have a diagonal letter
    if row + col == 6 and diagonal_letter:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid):
                return True
            grid[row][col] = ''
        return False

    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''

    return False

# Initial grid
grid = [
    ['c','','','','','',''],
    ['','','e','','','',''],
    ['d','e','g','','f','c','b'],
    ['e','g','','f','c','','d'],
    ['','a','f','','','',''],
    ['','f','c','b','','e','g'],
    ['f','c','','','e','g','']
]

# Verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            temp = grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, temp):
                valid = False
            grid[i][j] = temp

if valid and solve(grid):
    print_grid(grid)
else:
    print("No solution exists")