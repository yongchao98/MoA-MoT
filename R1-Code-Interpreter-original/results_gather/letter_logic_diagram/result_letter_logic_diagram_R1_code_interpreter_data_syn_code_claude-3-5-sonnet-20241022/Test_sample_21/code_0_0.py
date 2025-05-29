def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:  # Minor diagonal condition
        if letter != 'b':
            return False
            
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, only try 'b'
    if row + col == 6:
        if is_valid(grid, row, col, 'b'):
            grid[row][col] = 'b'
            if solve(grid):
                return True
            grid[row][col] = ''
        return False
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'd', 'g', 'f', 'a', 'b', ''],
    ['d', '', 'f', 'a', '', 'c', 'e'],
    ['g', 'f', '', 'b', 'c', 'e', ''],
    ['f', 'a', '', 'c', 'e', 'd', ''],
    ['a', 'b', '', '', 'd', 'g', 'f'],
    ['', 'c', '', 'd', '', '', 'a'],
    ['', 'e', 'd', 'g', 'f', 'a', 'b']
]

if solve(initial):
    for row in initial:
        print(','.join(row))
else:
    print("No solution exists")