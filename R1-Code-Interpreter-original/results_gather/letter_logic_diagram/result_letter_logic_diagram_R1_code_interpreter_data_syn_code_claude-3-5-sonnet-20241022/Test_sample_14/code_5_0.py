def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_minor_diagonal_letter(grid):
    # Get the letter that appears most on the minor diagonal in the initial grid
    letters = {}
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            letters[grid[i][j]] = letters.get(grid[i][j], 0) + 1
    return max(letters.items(), key=lambda x: x[1])[0] if letters else None

def verify_position(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve(grid, initial_grid, minor_diag_letter):
    # Find empty cell
    found = False
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                row, col = i, j
                found = True
                break
        if found:
            break
    
    if not found:
        return True

    # If this is a pre-filled position, skip it
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        if solve(grid, initial_grid, minor_diag_letter):
            return True
        grid[row][col] = ''
        return False

    # If this is on minor diagonal, must use minor_diag_letter
    if row + col == 6:
        if verify_position(grid, row, col, minor_diag_letter):
            grid[row][col] = minor_diag_letter
            if solve(grid, initial_grid, minor_diag_letter):
                return True
        grid[row][col] = ''
        return False

    # Try each letter
    for letter in 'abcdefg':
        if verify_position(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, minor_diag_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

grid = [row[:] for row in initial_grid]

# First, determine the letter that should be on the minor diagonal
minor_diag_letter = get_minor_diagonal_letter(initial_grid)
if not minor_diag_letter:
    minor_diag_letter = 'a'  # If no letter is found, start with 'a'

# Fill the minor diagonal first
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = minor_diag_letter

if solve(grid, initial_grid, minor_diag_letter):
    print_grid(grid)
else:
    print("No solution exists")