def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_fixed_diagonal_letter(grid):
    # Looking at pre-filled positions on diagonal
    diagonal_positions = [(i, 6-i) for i in range(7)]
    letters = {}
    for i, j in diagonal_positions:
        if grid[i][j] != '':
            letters[grid[i][j]] = letters.get(grid[i][j], 0) + 1
    # Return the letter that appears most on diagonal
    if letters:
        return max(letters.items(), key=lambda x: x[1])[0]
    return None

def get_available_letters(grid, row, col, diagonal_letter):
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    
    # If this is a diagonal position, only diagonal_letter is allowed
    if row + col == 6:
        return {diagonal_letter} if diagonal_letter not in (used_in_row | used_in_col) else set()
    
    # Otherwise, any unused letter except diagonal_letter (if this would complete a row/col)
    available = set('abcdefg') - used_in_row - used_in_col
    return available

def find_most_constrained_cell(grid, diagonal_letter):
    min_options = float('inf')
    best_cell = None
    best_options = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                options = get_available_letters(grid, i, j, diagonal_letter)
                if len(options) < min_options:
                    min_options = len(options)
                    best_cell = (i, j)
                    best_options = options
                    if min_options == 1:  # Can't get better than this
                        return best_cell, best_options
    
    return best_cell, best_options

def solve(grid, diagonal_letter):
    # Check if grid is complete
    if all(all(cell != '' for cell in row) for row in grid):
        return True
    
    # Find cell with fewest options
    cell, options = find_most_constrained_cell(grid, diagonal_letter)
    if not cell or not options:
        return False
    
    row, col = cell
    for letter in options:
        grid[row][col] = letter
        if solve(grid, diagonal_letter):
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

# First, determine the diagonal letter (must be 'e' based on initial grid)
diagonal_letter = 'e'

# Verify diagonal letter is compatible with pre-filled positions
diagonal_valid = True
for i in range(7):
    if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != diagonal_letter:
        diagonal_valid = False
        break

if diagonal_valid:
    # Fill in all diagonal positions with 'e'
    for i in range(7):
        if initial_grid[i][6-i] == '':
            initial_grid[i][6-i] = diagonal_letter
    
    if solve(initial_grid, diagonal_letter):
        print_grid(initial_grid)
    else:
        print("No solution found")
else:
    print("Invalid diagonal letter")