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
    if row + col == 6:  # Minor diagonal position
        # Check if it matches other filled cells on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
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
    
    # If it's on minor diagonal, we need to match existing diagonal letters
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', 'e', 'b', '', '', 'f', ''],
    ['', '', 'd', 'c', 'f', 'g', 'a'],
    ['', '', 'c', '', '', '', 'e'],
    ['', 'c', 'f', 'g', 'a', 'e', ''],
    ['c', '', 'g', '', '', '', 'd'],
    ['f', 'g', 'a', 'e', 'b', 'd', 'c'],
    ['', '', 'e', 'b', '', 'c', 'f']
]

if solve(initial_grid):
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")