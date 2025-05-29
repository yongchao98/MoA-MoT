def verify_grid(grid, initial):
    # Check initial constraints
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # Row check
            return False
        col = [grid[r][i] for r in range(7)]
        if len(set(col)) != 7:  # Column check
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:  # All should be same
        return False
    
    return True

def get_minor_diagonal_letter(initial):
    # Find the first letter on minor diagonal from initial grid
    for i in range(7):
        if initial[i][6-i]:
            return initial[i][6-i]
    return 'a'  # Default if none found

def get_valid_letters(grid, row, col, minor_diag_letter):
    if row + col == 6:  # If on minor diagonal
        return {minor_diag_letter}
    
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    return set('abcdefg') - used_row - used_col

def solve(grid, initial, minor_diag_letter, pos=0):
    if pos == 49:
        return verify_grid(grid, initial)
    
    row = pos // 7
    col = pos % 7
    
    # If this is a pre-filled position, skip to next
    if initial[row][col]:
        grid[row][col] = initial[row][col]
        return solve(grid, initial, minor_diag_letter, pos + 1)
    
    # Get valid letters for this position
    valid_letters = get_valid_letters(grid, row, col, minor_diag_letter)
    
    # Try each valid letter
    for letter in sorted(valid_letters):  # Sort for deterministic behavior
        grid[row][col] = letter
        if solve(grid, initial, minor_diag_letter, pos + 1):
            return True
        grid[row][col] = ''
    
    return False

# Initialize the puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# First, determine the minor diagonal letter
minor_diag_letter = get_minor_diagonal_letter(initial)

# Fill in initial values
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

# Fill minor diagonal with the determined letter
for i in range(7):
    if not initial[i][6-i]:
        grid[i][6-i] = minor_diag_letter

if solve(grid, initial, minor_diag_letter):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")