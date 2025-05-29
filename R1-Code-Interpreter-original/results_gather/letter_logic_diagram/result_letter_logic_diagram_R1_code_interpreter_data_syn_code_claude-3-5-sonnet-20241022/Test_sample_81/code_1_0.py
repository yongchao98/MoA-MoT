def verify_initial(initial):
    # First, collect all pre-filled positions
    prefilled = []
    diagonal_values = set()
    
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                prefilled.append((i, j, initial[i][j]))
                if i + j == 6:  # on minor diagonal
                    diagonal_values.add(initial[i][j])
    
    # Check if diagonal values are consistent
    if len(diagonal_values) > 1:
        return False
    
    # Check if pre-filled positions violate row/column constraints
    for i, j, val in prefilled:
        # Check row
        row_vals = [initial[i][k] for k in range(7) if k != j and initial[i][k] != '']
        if val in row_vals:
            return False
        
        # Check column
        col_vals = [initial[k][j] for k in range(7) if k != i and initial[k][j] != '']
        if val in col_vals:
            return False
    
    return True

def is_valid(grid, row, col, letter, initial):
    # Must match pre-filled value
    if initial[row][col] != '' and initial[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
            if initial[i][j] != '' and initial[i][j] != letter:
                return False
    
    return True

def solve(grid, initial, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If pre-filled, use that value and move to next position
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, pos + 1)
    
    # Get diagonal value if exists
    diagonal_letter = None
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if initial[i][j] != '':
                diagonal_letter = initial[i][j]
                break
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
    
    # Try letters
    letters = [diagonal_letter] if diagonal_letter else 'abcdefg'
    for letter in letters:
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'b', 'd', '', 'e', 'g', 'f'],
    ['', '', 'a', 'e', 'g', '', ''],
    ['d', 'a', '', '', '', 'c', ''],
    ['a', 'e', '', '', '', 'b', ''],
    ['', '', 'f', 'c', 'b', '', 'a'],
    ['', 'f', 'c', '', 'd', '', ''],
    ['f', '', 'b', '', '', '', 'g']
]

# First verify initial configuration is valid
print('<<<')
if not verify_initial(initial):
    print("Initial configuration is invalid!")
else:
    # Create empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    if solve(grid, initial):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")
print('>>>')