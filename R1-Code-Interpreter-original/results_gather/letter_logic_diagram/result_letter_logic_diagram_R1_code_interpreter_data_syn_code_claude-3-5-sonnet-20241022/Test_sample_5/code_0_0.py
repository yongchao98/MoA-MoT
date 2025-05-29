def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # If cell is on minor diagonal
        first_diag = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_diag = grid[i][6-i]
                break
        if first_diag and letter != first_diag:
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
    
    # If on minor diagonal, only try the letter that should be there
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'f', 'b', 'e', 'g', 'd', 'a'],
    ['', 'b', 'e', 'g', '', 'a', 'c'],
    ['', '', '', 'd', 'a', '', ''],
    ['e', '', '', 'a', 'c', 'f', ''],
    ['g', '', 'a', 'c', 'f', 'b', ''],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', '', '', 'e', 'g', 'd']
]

if solve_grid(initial):
    result = []
    for row in initial:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")