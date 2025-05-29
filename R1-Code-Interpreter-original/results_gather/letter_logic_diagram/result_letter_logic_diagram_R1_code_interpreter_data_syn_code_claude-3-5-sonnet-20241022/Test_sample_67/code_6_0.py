def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal - must be 'f'
    if row + col == 6 and letter != 'f':
        return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    # First fill all diagonal positions with 'f'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'f'
    
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', 'f'],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', 'f', '', 'c'],
    ['', '', '', 'f', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

if solve_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")