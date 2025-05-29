def is_valid_row(grid, row):
    letters = [x for x in grid[row] if x]
    return len(letters) == len(set(letters))

def is_valid_col(grid, col):
    letters = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(letters) == len(set(letters))

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If position is pre-filled, move to next position
    if grid[row][col]:
        return solve(grid, row, col + 1)
    
    # If on minor diagonal, must place 'g'
    if row + col == 6:
        grid[row][col] = 'g'
        if is_valid_row(grid, row) and is_valid_col(grid, col):
            if solve(grid, row, col + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        # Skip 'g' if this would exceed the number of 'g's allowed in row/column
        if letter == 'g' and row + col != 6:
            continue
            
        grid[row][col] = letter
        if is_valid_row(grid, row) and is_valid_col(grid, col):
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