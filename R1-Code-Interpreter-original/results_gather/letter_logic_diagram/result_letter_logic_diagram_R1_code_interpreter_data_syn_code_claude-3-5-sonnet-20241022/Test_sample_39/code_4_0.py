def is_valid(grid, row, col, letter, initial):
    # MUST match pre-filled cell
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, initial, pos=0):
    if pos == 49:
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If pre-filled, must use that value
    if initial[row][col] != '':
        if pos == 0 or is_valid(grid, row, col, initial[row][col], initial):
            grid[row][col] = initial[row][col]
            if solve(grid, initial, pos + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with strict pre-filled values
initial = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')