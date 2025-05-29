def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        first_diag = None
        for i in range(7):
            if grid[i][6-i] and grid[i][6-i] != letter:
                return False
            elif grid[i][6-i]:
                first_diag = grid[i][6-i]
    
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
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['', '', 'c', 'f', '', '', ''],
    ['e', 'c', 'f', '', '', '', ''],
    ['', 'f', 'g', 'b', 'a', '', 'e'],
    ['f', 'g', 'b', 'a', '', '', 'c'],
    ['g', 'b', 'a', 'd', 'e', '', 'f'],
    ['', '', '', '', 'c', '', ''],
    ['a', '', 'e', '', '', '', '']
]

if solve(initial):
    result = []
    for row in initial:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")