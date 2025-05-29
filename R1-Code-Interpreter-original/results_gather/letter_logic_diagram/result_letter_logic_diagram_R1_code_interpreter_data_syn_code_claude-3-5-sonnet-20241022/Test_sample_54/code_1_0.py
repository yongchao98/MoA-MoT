def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be 'd'
    if row + col == 6:  # Position is on minor diagonal
        if letter != 'd':
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
    
    # If position is on minor diagonal, only try 'd'
    if row + col == 6:
        letters = 'd'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', '', '', '', 'd', ''],
    ['', 'g', '', '', '', 'd', 'f'],
    ['', '', '', '', 'f', 'e', 'a'],
    ['b', '', '', '', 'f', 'e', ''],
    ['c', '', '', '', 'a', 'g', ''],
    ['d', 'f', '', '', 'g', '', ''],
    ['', 'e', 'a', '', 'b', '', '']
]

# First, ensure all minor diagonal positions are 'd'
for i in range(7):
    j = 6 - i
    if initial_grid[i][j] != '' and initial_grid[i][j] != 'd':
        print("Invalid input: minor diagonal must be 'd'")
        exit()
    initial_grid[i][j] = 'd'

# Solve the puzzle
if solve(initial_grid):
    # Format output
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")