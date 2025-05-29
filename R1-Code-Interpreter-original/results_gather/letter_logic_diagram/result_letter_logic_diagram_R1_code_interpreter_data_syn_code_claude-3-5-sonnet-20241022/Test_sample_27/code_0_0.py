def is_valid_row(grid, row):
    # Check if row contains each letter exactly once
    letters = [x for x in grid[row] if x]
    return len(letters) == len(set(letters))

def is_valid_col(grid, col):
    # Check if column contains each letter exactly once
    letters = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(letters) == len(set(letters))

def is_valid_minor_diag(grid):
    # Check if minor diagonal has same letter
    diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    return len(set(diag)) <= 1

def is_valid_placement(grid, row, col, letter):
    # Try placing letter and check constraints
    grid[row][col] = letter
    valid = is_valid_row(grid, row) and is_valid_col(grid, col)
    if row + col == 6:  # If on minor diagonal
        valid = valid and is_valid_minor_diag(grid)
    grid[row][col] = ''
    return valid

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    if grid[row][col]:  # Skip pre-filled cells
        return solve(grid, row, col + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
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

if solve(initial):
    result = []
    for row in initial:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution found")