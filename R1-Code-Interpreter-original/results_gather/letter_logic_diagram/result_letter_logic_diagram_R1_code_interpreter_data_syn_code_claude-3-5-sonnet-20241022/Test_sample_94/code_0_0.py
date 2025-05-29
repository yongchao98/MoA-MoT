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
        # Get the first non-empty letter on minor diagonal
        first_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_letter = grid[i][6-i]
                break
        if first_letter and letter != first_letter:
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
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['a','f','g','c','e','d','b'],
    ['f','','c','','','b','a'],
    ['g','c','e','','','a','f'],
    ['c','','','','','f',''],
    ['e','d','','a','f','',''],
    ['d','b','a','','g','','e'],
    ['b','','','','','','d']
]

# Solve the puzzle
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