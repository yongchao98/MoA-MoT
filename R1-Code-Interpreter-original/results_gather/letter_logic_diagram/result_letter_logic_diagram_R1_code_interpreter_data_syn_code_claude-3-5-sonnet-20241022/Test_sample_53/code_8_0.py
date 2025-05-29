def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'c'
    if row + col == 6 and letter != 'c':
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
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If on minor diagonal, only try 'c'
    letters = ['c'] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'f', '', '', '', 'c'],  # Force 'c' on minor diagonal
    ['', '', 'd', 'e', 'g', 'c', ''],
    ['', '', 'e', '', 'c', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', 'c', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['c', '', '', 'f', 'd', 'e', 'g']
]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")