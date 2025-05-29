def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_possible_diagonal_letters(grid):
    # Get all letters that appear on the diagonal in the initial grid
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    return diagonal_letters

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    # Check diagonal
    if row + col == 6 and letter != diag_letter:
        return False
    
    return True

def solve(grid, diag_letter, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        # Verify diagonal constraint for pre-filled cells
        if row + col == 6 and grid[row][col] != diag_letter:
            return False
        return solve(grid, diag_letter, row, col + 1)
    
    # For diagonal positions, only try the diagonal letter
    if row + col == 6:
        if is_valid(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal positions, try all letters
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diag_letter):
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

# Try each letter as the diagonal letter
for diag_letter in 'abcdefg':
    # Create a copy of the initial grid
    grid_copy = [row[:] for row in initial_grid]
    
    # Check if this diagonal letter is compatible with pre-filled positions
    valid = True
    for i in range(7):
        if grid_copy[i][6-i] != '' and grid_copy[i][6-i] != diag_letter:
            valid = False
            break
    
    if valid and solve(grid_copy, diag_letter):
        print_grid(grid_copy)
        break