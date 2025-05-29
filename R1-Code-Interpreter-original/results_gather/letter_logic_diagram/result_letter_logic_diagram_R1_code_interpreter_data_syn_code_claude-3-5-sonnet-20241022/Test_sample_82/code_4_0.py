def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_valid_position(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    # Check column        
    for x in range(7):
        if grid[x][col] == num:
            return False
    return True

def solve(grid, diagonal_letter, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True

    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row, col + 1)

    # If this is a diagonal position, we must use diagonal_letter
    if row + col == 6:
        if check_valid_position(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
        grid[row][col] = ''
        return False

    # Try each possible letter for non-diagonal positions
    for letter in 'abcdefg':
        if check_valid_position(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['f','','e','','','',''],
    ['','e','','d','','','f'],
    ['','','','','g','f','a'],
    ['','d','b','','','a',''],
    ['d','b','','','a','e',''],
    ['b','g','','a','e','c','d'],
    ['','f','a','','','d','']
]

# First, check which letters appear on the diagonal
diagonal_cells = [(i, 6-i) for i in range(7)]
diagonal_letters = set()
for i, j in diagonal_cells:
    if initial_grid[i][j] != '':
        diagonal_letters.add(initial_grid[i][j])

# If we have a pre-filled diagonal letter, use it; otherwise try each possible letter
if len(diagonal_letters) > 0:
    diagonal_letter = list(diagonal_letters)[0]
else:
    diagonal_letter = 'c'  # Start with 'c' as it appears in one position

# Create a working copy of the grid
grid = [row[:] for row in initial_grid]

# Fill all diagonal positions with the chosen letter
for i, j in diagonal_cells:
    if grid[i][j] == '':
        grid[i][j] = diagonal_letter

# Try to solve the grid
if solve(grid, diagonal_letter):
    print_grid(grid)
else:
    # If failed, try other diagonal letters
    for new_diagonal_letter in 'abcdefg':
        if new_diagonal_letter == diagonal_letter:
            continue
        grid = [row[:] for row in initial_grid]
        # Fill diagonal with new letter
        for i, j in diagonal_cells:
            if grid[i][j] == '':
                grid[i][j] = new_diagonal_letter
        if solve(grid, new_diagonal_letter):
            print_grid(grid)
            break