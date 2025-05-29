def find_diagonal_letter(grid):
    # First check if any position on minor diagonal is filled
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid_with_diagonal(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if position is on minor diagonal
    if row + col == 6 and letter != diag_letter:
        return False
    
    return True

def solve_with_fixed_diagonal(grid, diag_letter):
    # First fill all minor diagonal positions
    grid_copy = [row[:] for row in grid]
    for i in range(7):
        if grid_copy[i][6-i] == '':
            grid_copy[i][6-i] = diag_letter
        elif grid_copy[i][6-i] != diag_letter:
            return False
    
    if solve_remaining(grid_copy, diag_letter):
        for i in range(7):
            grid[i] = grid_copy[i][:]
        return True
    return False

def solve_remaining(grid, diag_letter):
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
        if is_valid_with_diagonal(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve_remaining(grid, diag_letter):
                return True
            grid[row][col] = ''
    
    return False

def solve_grid(grid):
    # First find if any diagonal position is filled
    diag_letter = find_diagonal_letter(grid)
    
    if diag_letter:
        return solve_with_fixed_diagonal(grid, diag_letter)
    else:
        # Try each possible letter for diagonal
        for letter in 'abcdefg':
            if solve_with_fixed_diagonal(grid, letter):
                return True
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