def get_valid_letters(grid, row, col, initial):
    # If pre-filled, only that letter is valid
    if initial[row][col] != '':
        return {initial[row][col]}
    
    # Start with all possible letters
    valid = set('abcdefg')
    
    # Remove letters from same row
    for j in range(7):
        if grid[row][j] != '':
            valid.discard(grid[row][j])
        # Also check pre-filled positions
        if initial[row][j] != '' and j != col:
            valid.discard(initial[row][j])
    
    # Remove letters from same column
    for i in range(7):
        if grid[i][col] != '':
            valid.discard(grid[i][col])
        # Also check pre-filled positions
        if initial[i][col] != '' and i != row:
            valid.discard(initial[i][col])
    
    # If on minor diagonal, must match existing diagonal letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
            if initial[i][j] != '':
                diagonal_letter = initial[i][j]
                break
        if diagonal_letter:
            valid = {diagonal_letter} & valid
    
    return valid

def verify_solution(grid, initial):
    # Verify pre-filled positions
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    
    # Verify rows and columns
    for i in range(7):
        row = set(grid[i])
        col = set(grid[r][i] for r in range(7))
        if len(row) != 7 or len(col) != 7:
            return False
    
    # Verify minor diagonal
    diagonal = set(grid[i][6-i] for i in range(7))
    if len(diagonal) != 1:
        return False
    
    return True

def solve(grid, initial):
    # Find empty cell with minimum valid choices
    min_choices = 8
    min_pos = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                valid = get_valid_letters(grid, i, j, initial)
                if len(valid) == 0:  # No valid choices
                    return False
                if len(valid) < min_choices:
                    min_choices = len(valid)
                    min_pos = (i, j, valid)
    
    if not min_pos:  # Grid is filled
        return verify_solution(grid, initial)
    
    row, col, valid_letters = min_pos
    
    # Try each valid letter
    for letter in valid_letters:
        grid[row][col] = letter
        if solve(grid, initial):
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

# Create working grid and copy initial values
grid = [['' for _ in range(7)] for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j] != '':
            grid[i][j] = initial[i][j]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')