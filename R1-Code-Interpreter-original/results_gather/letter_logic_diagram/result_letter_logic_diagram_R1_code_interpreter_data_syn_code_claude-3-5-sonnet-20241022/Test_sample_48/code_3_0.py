def is_valid_move(grid, row, col, letter):
    # If on minor diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal cells
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If on minor diagonal, only try 'g'
    if row + col == 6:
        letters = ['g']
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', ''],
    ['', '', 'c', 'g', 'f', '', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', '', 'a'],
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")