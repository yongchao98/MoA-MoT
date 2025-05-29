def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_prefilled_diagonal(initial):
    # Get all pre-filled letters on the minor diagonal
    diag_letters = []
    for i in range(7):
        if initial[i][6-i]:
            diag_letters.append((i, 6-i, initial[i][6-i]))
    return diag_letters

def is_valid_move(grid, row, col, letter, initial):
    # Must match initial configuration
    if initial[row][col] and initial[row][col] != letter:
        return False
    
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal
    if row + col == 6:
        diag = [(i, 6-i) for i in range(7)]
        for i, j in diag:
            if grid[i][j] and grid[i][j] != letter:
                return False
    
    return True

def copy_initial_state(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j]:
                grid[i][j] = initial[i][j]

def solve(grid, initial, pos=0):
    if pos >= 49:  # 7x7 = 49 cells
        return True
    
    row, col = pos // 7, pos % 7
    
    # Skip pre-filled cells
    if initial[row][col]:
        return solve(grid, initial, pos + 1)
    
    # If on minor diagonal, must match other diagonal letters
    if row + col == 6:
        diag_letters = set()
        for i in range(7):
            if grid[i][6-i]:
                diag_letters.add(grid[i][6-i])
        if diag_letters:
            letters = diag_letters
        else:
            letters = set('abcdefg')
    else:
        letters = set('abcdefg')
    
    # Try each possible letter
    for letter in letters:
        if is_valid_move(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy initial state
copy_initial_state(grid, initial)

# Get pre-filled diagonal positions
diag_positions = get_prefilled_diagonal(initial)

# Try each possible letter for the diagonal
for diag_letter in 'abcdefg':
    test_grid = [row[:] for row in grid]
    valid = True
    
    # Check if this diagonal letter is compatible with pre-filled diagonal positions
    for _, _, letter in diag_positions:
        if letter != diag_letter:
            valid = False
            break
    
    if not valid:
        continue
        
    # Fill all diagonal positions with this letter
    for i in range(7):
        if test_grid[i][6-i] == '':
            if is_valid_move(test_grid, i, 6-i, diag_letter, initial):
                test_grid[i][6-i] = diag_letter
            else:
                valid = False
                break
    
    if valid and solve(test_grid, initial):
        print_solution(test_grid)
        break