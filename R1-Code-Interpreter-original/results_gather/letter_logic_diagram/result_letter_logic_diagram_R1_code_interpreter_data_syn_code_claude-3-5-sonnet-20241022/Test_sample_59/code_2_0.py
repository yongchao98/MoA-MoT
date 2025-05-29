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
    minor_diag_letter = get_minor_diagonal_letter(grid)
    if row + col == 6:  # If position is on minor diagonal
        if minor_diag_letter and letter != minor_diag_letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal positions
    minor_diag_letter = get_minor_diagonal_letter(grid)
    if minor_diag_letter:
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    # Get minor diagonal letter
    minor_diag_letter = get_minor_diagonal_letter(grid)
    if not minor_diag_letter:
        # If no letter is set for minor diagonal, we need to try each possibility
        for letter in 'abcdefg':
            grid_copy = [row[:] for row in grid]
            # Fill first empty position in minor diagonal
            for i in range(7):
                if grid_copy[i][6-i] == '':
                    grid_copy[i][6-i] = letter
                    break
            if solve_with_fixed_diagonal(grid_copy):
                for i in range(7):
                    grid[i] = grid_copy[i]
                return True
        return False
    else:
        return solve_with_fixed_diagonal(grid)

def solve_with_fixed_diagonal(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    if row + col == 6:  # If on minor diagonal
        minor_letter = get_minor_diagonal_letter(grid)
        if is_valid(grid, row, col, minor_letter):
            grid[row][col] = minor_letter
            if solve_with_fixed_diagonal(grid):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_with_fixed_diagonal(grid):
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