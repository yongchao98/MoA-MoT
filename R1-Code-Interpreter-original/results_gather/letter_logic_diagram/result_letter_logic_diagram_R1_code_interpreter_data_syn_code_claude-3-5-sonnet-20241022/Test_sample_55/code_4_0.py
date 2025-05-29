def fill_diagonal_and_check(grid):
    # Fill all minor diagonal cells with 'b' and check if it conflicts with pre-filled cells
    for i in range(7):
        j = 6-i
        if grid[i][j] and grid[i][j] != 'b':
            return False
        grid[i][j] = 'b'
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    return True

def find_empty(grid):
    # Fill non-diagonal cells row by row
    for i in range(7):
        for j in range(7):
            if j != 6-i and not grid[i][j]:  # Skip diagonal cells
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Try letters in different order
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['a', 'g', '', '', 'f', '', 'b'],
    ['', 'd', '', '', 'e', '', 'a'],
    ['d', '', 'f', 'e', '', 'a', 'g'],
    ['c', 'f', 'e', 'b', 'a', 'g', ''],
    ['', '', 'b', '', 'g', 'd', 'c'],
    ['e', 'b', 'a', 'g', 'd', 'c', ''],
    ['b', '', 'g', 'd', 'c', 'f', 'e']
]

# Make a copy of the original grid
original_grid = [row[:] for row in grid]

# First fill all diagonal cells with 'b'
if fill_diagonal_and_check(grid):
    if solve(grid):
        result = []
        for row in grid:
            result.append(','.join(row))
        print('<<<')
        for row in result:
            print(row)
        print('>>>')
    else:
        print("No solution exists")
else:
    print("Invalid diagonal configuration")