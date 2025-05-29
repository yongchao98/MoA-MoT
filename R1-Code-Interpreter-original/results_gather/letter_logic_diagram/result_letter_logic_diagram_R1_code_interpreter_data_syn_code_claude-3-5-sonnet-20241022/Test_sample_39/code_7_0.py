def initialize_grid(initial):
    grid = [['' for _ in range(7)] for _ in range(7)]
    # Copy initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
    return grid

def is_valid_state(grid, row, col, letter, initial):
    # Must match initial configuration
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, initial):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If this position is on the minor diagonal, we need to check existing diagonal values
    diagonal_letter = None
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
    
    # Try each possible letter
    letters = [diagonal_letter] if diagonal_letter else 'abcdefg'
    for letter in letters:
        if is_valid_state(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initial configuration
initial = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

# Create and solve grid
grid = initialize_grid(initial)

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')