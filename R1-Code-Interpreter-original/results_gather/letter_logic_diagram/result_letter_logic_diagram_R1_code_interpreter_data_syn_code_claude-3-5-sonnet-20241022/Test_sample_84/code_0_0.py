def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:  # Cell is on minor diagonal
        # Find any existing letter on minor diagonal
        for i in range(7):
            if grid[i][6-i] and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, use the existing diagonal letter if any exists
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'f', '', '', 'e', 'c'],
    ['d', 'f', 'b', 'g', 'e', 'c', ''],
    ['f', 'b', 'g', 'e', '', 'a', ''],
    ['', 'g', '', '', 'a', 'd', ''],
    ['', 'e', '', '', 'd', '', 'b'],
    ['e', '', 'a', 'd', 'f', 'b', 'g'],
    ['', '', 'd', 'f', '', 'g', '']
]

if solve_grid(initial_grid):
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")