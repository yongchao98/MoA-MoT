def get_minor_diagonal_letter(grid):
    # Get the letter that must be on minor diagonal from pre-filled positions
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # If position is on minor diagonal
        diag_letter = get_minor_diagonal_letter(grid)
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid):
    # First, determine and fill the minor diagonal
    diag_letter = get_minor_diagonal_letter(grid)
    if not diag_letter:
        # Try each possible letter for the diagonal
        for letter in 'abcdefg':
            grid_copy = [row[:] for row in grid]
            # Fill entire diagonal with this letter
            valid = True
            for i in range(7):
                if grid_copy[i][6-i] != '' and grid_copy[i][6-i] != letter:
                    valid = False
                    break
                grid_copy[i][6-i] = letter
            
            if valid and solve_remaining(grid_copy):
                for i in range(7):
                    grid[i] = grid_copy[i]
                return True
        return False
    else:
        # Fill entire diagonal with the known letter
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diag_letter
        return solve_remaining(grid)

def solve_remaining(grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")