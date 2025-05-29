def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    return True

def get_diagonal_options(grid):
    # Get all possible letters that could be on the diagonal
    # based on pre-filled diagonal positions
    options = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            options = {grid[i][6-i]}
            break
    return options

def solve(grid, diag_letter, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, diag_letter, row, col + 1)
    
    letters = [diag_letter] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, diag_letter, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', '', 'a', ''],
    ['', 'e', '', 'd', '', '', ''],
    ['', 'c', '', 'a', '', 'g', 'f'],
    ['', '', '', 'b', '', 'f', 'e'],
    ['d', 'a', '', '', '', 'e', ''],
    ['', '', '', 'f', 'e', '', ''],
    ['b', 'g', 'f', 'e', '', '', '']
]

# Try each possible diagonal letter
for diag_letter in get_diagonal_options(initial_grid):
    grid_copy = [row[:] for row in initial_grid]
    # First, fill all diagonal positions with the chosen letter
    for i in range(7):
        if grid_copy[i][6-i] == '':
            if is_valid(grid_copy, i, 6-i, diag_letter):
                grid_copy[i][6-i] = diag_letter
            else:
                break
    else:
        if solve(grid_copy, diag_letter):
            print_grid(grid_copy)
            break