def print_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_row(grid, row):
    letters = [x for x in grid[row] if x]
    return len(letters) == len(set(letters))

def is_valid_col(grid, col):
    letters = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(letters) == len(set(letters))

def check_initial_grid(grid):
    # Check minor diagonal for consistency
    diag_letter = None
    for i in range(7):
        if grid[i][6-i]:
            if diag_letter is None:
                diag_letter = grid[i][6-i]
            elif grid[i][6-i] != diag_letter:
                return False, None
    return True, diag_letter

def solve(grid, diag_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, diag_letter, row + 1, 0)
    
    # If position is pre-filled, move to next position
    if grid[row][col]:
        return solve(grid, diag_letter, row, col + 1)
    
    # If on minor diagonal, must use diagonal letter
    if row + col == 6:
        grid[row][col] = diag_letter
        if is_valid_row(grid, row) and is_valid_col(grid, col):
            if solve(grid, diag_letter, row, col + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each letter except diagonal letter
    available_letters = [l for l in 'abcdefg' if l != diag_letter]
    for letter in available_letters:
        grid[row][col] = letter
        if is_valid_row(grid, row) and is_valid_col(grid, col):
            if solve(grid, diag_letter, row, col + 1):
                return True
        grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'd', 'f', 'e', 'b', 'a', 'g'],
    ['', 'f', '', 'b', 'a', 'g', ''],
    ['', 'e', 'b', '', '', 'c', 'd'],
    ['', 'b', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

# First check if initial grid is valid and get diagonal letter
valid, diag_letter = check_initial_grid(initial)

if not valid:
    print("Invalid initial grid")
else:
    if solve(initial, diag_letter):
        print_grid(initial)
    else:
        print("No solution found")