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
    # Then other positions
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
    letters = 'c' if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

# First fill all minor diagonal positions with 'c'
for i in range(7):
    if initial_grid[i][6-i] == '':
        initial_grid[i][6-i] = 'c'

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")