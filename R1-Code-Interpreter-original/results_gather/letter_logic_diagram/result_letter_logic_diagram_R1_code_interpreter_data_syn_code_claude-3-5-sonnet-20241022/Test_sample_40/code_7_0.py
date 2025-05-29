def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_config():
    initial_grid = [
        ['', 'e', '', 'g', '', 'c', 'b'],
        ['', 'd', '', '', '', 'b', 'a'],
        ['', '', '', 'c', '', '', ''],
        ['', '', 'c', 'b', '', 'e', 'd'],
        ['f', '', '', '', 'e', 'd', ''],
        ['', '', '', '', '', '', ''],
        ['b', 'a', '', 'd', 'g', '', '']
    ]
    
    # Check for conflicts in initial values
    for i in range(7):
        row_letters = [x for x in initial_grid[i] if x != '']
        if len(row_letters) != len(set(row_letters)):
            return None
        
        col_letters = [initial_grid[j][i] for j in range(7) if initial_grid[j][i] != '']
        if len(col_letters) != len(set(col_letters)):
            return None
    
    # Find diagonal letter
    diag_letter = None
    diag_positions = []
    for i in range(7):
        if initial_grid[i][6-i] != '':
            if diag_letter is None:
                diag_letter = initial_grid[i][6-i]
            elif initial_grid[i][6-i] != diag_letter:
                return None
        diag_positions.append((i, 6-i))
    
    return initial_grid, diag_letter, diag_positions

def is_valid(grid, row, col, letter, diag_positions, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check diagonal
    if (row, col) in diag_positions and letter != diag_letter:
        return False
    
    return True

def solve(grid, pos, diag_positions, diag_letter):
    if pos == 49:
        return True
    
    row = pos // 7
    col = pos % 7
    
    # Skip if cell is pre-filled
    if grid[row][col] != '':
        return solve(grid, pos + 1, diag_positions, diag_letter)
    
    # If on diagonal, must use diagonal letter
    if (row, col) in diag_positions:
        if is_valid(grid, row, col, diag_letter, diag_positions, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, pos + 1, diag_positions, diag_letter):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diag_positions, diag_letter):
            grid[row][col] = letter
            if solve(grid, pos + 1, diag_positions, diag_letter):
                return True
            grid[row][col] = ''
    
    return False

# Check initial configuration and get necessary information
config = check_initial_config()
if config is None:
    print('<<<')
    print("No solution exists - initial configuration is invalid")
    print('>>>')
else:
    initial_grid, diag_letter, diag_positions = config
    
    # Create working grid
    grid = [row[:] for row in initial_grid]
    
    print('<<<')
    if solve(grid, 0, diag_positions, diag_letter):
        print_grid(grid)
    else:
        print("No solution exists")
    print('>>>')