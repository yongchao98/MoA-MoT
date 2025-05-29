def get_valid_letters(grid, row, col, initial):
    if initial[row][col] != '':
        return {initial[row][col]}
    
    # All possible letters minus those in row and column
    valid = set('abcdefg')
    
    # Remove letters from same row
    for j in range(7):
        if grid[row][j] != '':
            valid.discard(grid[row][j])
    
    # Remove letters from same column
    for i in range(7):
        if grid[i][col] != '':
            valid.discard(grid[i][col])
    
    # If on minor diagonal, must match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                return {grid[i][j]} & valid
    
    return valid

def solve(grid, initial):
    # Find cell with minimum possible valid letters
    min_choices = 8
    min_pos = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                valid = get_valid_letters(grid, i, j, initial)
                if len(valid) < min_choices:
                    min_choices = len(valid)
                    min_pos = (i, j, valid)
                    if min_choices == 1:  # Can't get better than this
                        break
        if min_choices == 1:
            break
    
    if not min_pos:  # Grid is filled
        return True
    
    row, col, valid_letters = min_pos
    
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
        grid[i][j] = initial[i][j]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')