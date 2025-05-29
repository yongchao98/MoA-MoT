def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_placement(grid, row, col, letter, initial):
    # Check initial constraint
    if initial[row][col] and initial[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return (i, j)
    return None

def solve_grid(grid, initial):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    available = set('abcdefg') - used_row - used_col
    
    for letter in sorted(available):  # Try letters in alphabetical order
        if is_valid_placement(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve_grid(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initialize puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Step 1: Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

# Step 2: Fill minor diagonal with 'a' (checking if it's possible)
can_use_a = True
for i in range(7):
    j = 6 - i
    if initial[i][j] and initial[i][j] != 'a':
        can_use_a = False
        break

if can_use_a:
    # Fill minor diagonal with 'a'
    for i in range(7):
        j = 6 - i
        if not grid[i][j]:
            grid[i][j] = 'a'

    # Step 3: Solve the rest of the grid
    if solve_grid(grid, initial):
        print_solution(grid)
    else:
        print("No solution exists")
else:
    print("Cannot use 'a' for minor diagonal")