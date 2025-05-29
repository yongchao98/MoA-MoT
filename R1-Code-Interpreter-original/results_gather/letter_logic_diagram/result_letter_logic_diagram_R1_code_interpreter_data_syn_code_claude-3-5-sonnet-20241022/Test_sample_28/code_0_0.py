def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal if this cell is on it
    if row + col == 6:
        diag = get_minor_diagonal(grid)
        filled_diag = [x for x in diag if x]
        if filled_diag and filled_diag[0] != num:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, must match other diagonal elements
    if row + col == 6:
        diag = get_minor_diagonal(grid)
        filled_diag = [x for x in diag if x]
        if filled_diag:
            letters = filled_diag[0]
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")