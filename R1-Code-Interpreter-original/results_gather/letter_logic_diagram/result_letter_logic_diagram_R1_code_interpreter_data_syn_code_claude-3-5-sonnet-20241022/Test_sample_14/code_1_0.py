def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_position(grid, row, col, letter, initial_grid):
    # Check if this violates initial configuration
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def get_minor_diagonal_letter(initial_grid):
    # Find the most frequent letter on the minor diagonal in initial grid
    counts = {}
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            counts[initial_grid[i][j]] = counts.get(initial_grid[i][j], 0) + 1
    
    if counts:
        return max(counts.items(), key=lambda x: x[1])[0]
    return 'b'  # Default to 'b' if no letter is found

def solve(grid, initial_grid, minor_diag_letter):
    def find_next_empty():
        # First fill minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j, True)
        
        # Then fill positions with initial values
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and initial_grid[i][j] != '':
                    return (i, j, False)
        
        # Then fill remaining positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j, False)
        return None

    pos = find_next_empty()
    if not pos:
        return True

    row, col, is_diagonal = pos

    # If this position has an initial value, try only that
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if verify_position(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, minor_diag_letter):
                return True
            grid[row][col] = ''
        return False

    # If this is on minor diagonal, use minor_diag_letter
    if is_diagonal:
        if verify_position(grid, row, col, minor_diag_letter, initial_grid):
            grid[row][col] = minor_diag_letter
            if solve(grid, initial_grid, minor_diag_letter):
                return True
            grid[row][col] = ''
        return False

    # Try all possible letters
    for letter in 'abcdefg':
        if verify_position(grid, row, col, letter, initial_grid):
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

# Create working grid
grid = [row[:] for row in initial_grid]

# First determine the letter for minor diagonal
minor_diag_letter = get_minor_diagonal_letter(initial_grid)

# Fill minor diagonal with the chosen letter
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = minor_diag_letter

if solve(grid, initial_grid, minor_diag_letter):
    print_grid(grid)
else:
    print("No solution exists")