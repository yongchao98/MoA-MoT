def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # if on minor diagonal
        return letter == 'g'
    elif letter == 'g':  # if not on minor diagonal, can't be 'g'
        return False
    
    return True

def find_empty(grid):
    # Find first empty cell, prioritizing minor diagonal positions
    # First check minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    
    # Then check other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and i + j != 6:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = ['g'] if row + col == 6 else 'abcdef'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
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

# Validate initial grid
valid = True
for i in range(7):
    if initial[i][6-i] and initial[i][6-i] != 'g':
        valid = False
        break

if not valid:
    print("Invalid initial grid")
else:
    if solve(initial):
        print('<<<')
        for row in initial:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution found")