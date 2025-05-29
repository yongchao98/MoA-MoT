def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal - must be 'f' (we'll choose 'f' as our diagonal letter)
    if row + col == 6 and letter != 'f':
        return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if i + (6-i) == 6 and not grid[i][6-i]:
            return i, 6-i
    # Then fill rest
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
    # If on minor diagonal, only try 'f'
    if row + col == 6:
        letters = ['f']
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Verify if 'f' is possible on diagonal by checking pre-filled cells
can_use_f = True
for i in range(7):
    if initial[i][6-i] and initial[i][6-i] != 'f':
        can_use_f = False
        break

if can_use_f:
    # Force 'f' on all diagonal positions
    for i in range(7):
        if not initial[i][6-i]:
            initial[i][6-i] = 'f'

# Solve
if solve(initial):
    # Print solution
    result = []
    for row in initial:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")