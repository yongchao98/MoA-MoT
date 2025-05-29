def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:
        # Get the letter that should be on the minor diagonal
        diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diag_letter = grid[i][j]
                break
        if diag_letter is not None and letter != diag_letter:
            return False
    
    return True

def find_empty(grid):
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
    letters = 'abcdefg'
    
    # If on minor diagonal, try the existing diagonal letter first
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j] + ''.join(l for l in 'abcdefg' if l != grid[i][j])
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
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

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")