def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_empty_cell(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        first_diag = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_diag = grid[i][6-i]
                break
        if first_diag and letter != first_diag:
            return False
    
    return True

def solve(grid):
    empty = get_empty_cell(grid)
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

# Initialize grid
grid = [
    ['d','e','c','a','f','','g'],
    ['','c','','','b','g',''],
    ['c','a','f','','','d','e'],
    ['a','','','g','d','e','c'],
    ['f','','g','d','e','c','a'],
    ['b','g','','','c','','f'],
    ['g','d','e','c','','f','']
]

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")