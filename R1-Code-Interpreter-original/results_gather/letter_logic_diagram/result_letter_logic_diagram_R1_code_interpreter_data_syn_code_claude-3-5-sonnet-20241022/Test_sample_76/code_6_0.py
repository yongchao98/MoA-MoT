def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_letter(initial):
    # Get the required letter for minor diagonal from initial grid
    for i in range(7):
        if initial[i][6-i] != '':
            return initial[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial):
    # Must match initial grid
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:
        diag_letter = get_diagonal_letter(initial)
        if diag_letter and letter != diag_letter:
            return False

    return True

def solve(grid, initial, pos=0):
    if pos >= 49:
        return True

    row = pos // 7
    col = pos % 7

    # If this is a pre-filled position, copy it and move on
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, pos + 1)

    # Get diagonal letter if we're on the diagonal
    diag_letter = None
    if row + col == 6:
        diag_letter = get_diagonal_letter(initial)

    # Try letters
    letters = [diag_letter] if diag_letter else list('abcdefg')
    for letter in letters:
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''

    return False

# Initialize with exact initial values
initial = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# First, verify the initial grid is valid
valid_initial = True
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            letter = initial[i][j]
            initial[i][j] = ''
            if not is_valid(initial, i, j, letter, initial):
                valid_initial = False
            initial[i][j] = letter

print('<<<')
if valid_initial and solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No valid solution exists")
print('>>>')