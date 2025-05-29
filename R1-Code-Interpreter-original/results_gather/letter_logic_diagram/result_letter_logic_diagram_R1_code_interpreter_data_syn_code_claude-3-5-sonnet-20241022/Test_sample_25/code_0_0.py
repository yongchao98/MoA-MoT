def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency if this cell is on minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
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
    
    # If on minor diagonal, only try the letter that's already there (if any)
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

# Initialize grid with given values
initial_grid = [
    ['f','','d','','','g','a'],
    ['','','e','','','',''],
    ['','e','','g','','','c'],
    ['','','','a','','c','d'],
    ['','g','','f','c','d','e'],
    ['','','','','','',''],
    ['','','c','d','e','','']
]

# Solve the puzzle
if solve_grid(initial_grid):
    # Print solution in required format
    result = ""
    for row in initial_grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")