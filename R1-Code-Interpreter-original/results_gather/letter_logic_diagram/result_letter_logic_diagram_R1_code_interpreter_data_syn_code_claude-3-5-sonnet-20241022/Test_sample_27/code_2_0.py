def is_valid_row(grid, row):
    letters = [x for x in grid[row] if x]
    return len(letters) == len(set(letters))

def is_valid_col(grid, col):
    letters = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(letters) == len(set(letters))

def is_valid_placement(grid, row, col, letter):
    # If on minor diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
        
    # Try placing and check row/column constraints
    grid[row][col] = letter
    valid = is_valid_row(grid, row) and is_valid_col(grid, col)
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
    
    # If on minor diagonal, only try 'g'
    if row + col == 6:
        if is_valid_placement(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

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