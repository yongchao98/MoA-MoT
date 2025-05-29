def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal - must be 'e'
    if row + col == 6 and letter != 'e':
        return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def get_valid_letters(grid, row, col, used_in_row, used_in_col):
    if row + col == 6:  # Minor diagonal must be 'e'
        return ['e'] if 'e' not in used_in_row and 'e' not in used_in_col else []
    
    valid = []
    for letter in 'abcdefg':
        if (letter not in used_in_row and 
            letter not in used_in_col and 
            is_valid(grid, row, col, letter)):
            valid.append(letter)
    return valid

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Get letters already used in this row and column
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    
    valid_letters = get_valid_letters(grid, row, col, used_in_row, used_in_col)
    
    for letter in valid_letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['e','c','b','','f','d','g'],
    ['c','','','f','','','e'],
    ['b','a','','','g','e','c'],
    ['','f','','','e','c',''],
    ['','d','','e','','',''],
    ['','','e','c','b','','f'],
    ['','','c','','','','']
]

# First, verify and fill minor diagonal with 'e'
for i in range(7):
    if grid[i][6-i] and grid[i][6-i] != 'e':
        print("Invalid initial grid: minor diagonal must be 'e'")
        exit()
    grid[i][6-i] = 'e'

# Solve puzzle
if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")