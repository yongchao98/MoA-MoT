def get_minor_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
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
    
    # If this is on minor diagonal
    if row + col == 6:
        minor_letter = get_minor_diagonal_letter(grid)
        if minor_letter and letter != minor_letter:
            return False
    
    return True

def solve_grid(grid):
    # First, determine minor diagonal letter if any exists
    minor_letter = get_minor_diagonal_letter(grid)
    if not minor_letter:
        # If no minor diagonal letter is set, we need to try each possibility
        for letter in 'abcdefg':
            # Try this letter for all minor diagonal positions
            if try_minor_diagonal(grid, letter):
                return True
        return False
    else:
        # Fill minor diagonal with the existing letter
        if not try_minor_diagonal(grid, minor_letter):
            return False
        # Then fill the rest
        return solve_remaining(grid)

def try_minor_diagonal(grid, letter):
    # First fill all minor diagonal positions with the letter
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        j = 6 - i
        if temp_grid[i][j] == '':
            if not is_valid(temp_grid, i, j, letter):
                return False
            temp_grid[i][j] = letter
    
    # Then solve the rest
    return solve_remaining(temp_grid) and update_grid(grid, temp_grid)

def solve_remaining(grid):
    empty = find_empty(grid)
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

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def update_grid(grid, temp_grid):
    for i in range(7):
        for j in range(7):
            grid[i][j] = temp_grid[i][j]
    return True

# Initial grid
initial_grid = [
    ['', 'd', '', '', 'b', 'e', 'g'],
    ['', '', '', 'b', '', 'g', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution exists")