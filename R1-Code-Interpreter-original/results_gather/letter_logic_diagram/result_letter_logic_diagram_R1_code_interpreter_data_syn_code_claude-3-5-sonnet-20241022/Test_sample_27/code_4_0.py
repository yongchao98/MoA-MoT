def is_valid_row(grid, row):
    letters = [x for x in grid[row] if x]
    return len(letters) == len(set(letters))

def is_valid_col(grid, col):
    letters = [grid[r][col] for r in range(7) if grid[r][col]]
    return len(letters) == len(set(letters))

def fill_minor_diagonal(grid):
    # Fill all minor diagonal positions with 'g'
    for i in range(7):
        if not grid[i][6-i]:  # if empty
            grid[i][6-i] = 'g'
    return grid

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip if position is pre-filled or on minor diagonal
    if grid[row][col]:
        return solve(grid, row, col + 1)
    
    # Try each letter except 'g' (since 'g' is only for minor diagonal)
    for letter in 'abcdef':
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

# First fill all minor diagonal positions with 'g'
initial = fill_minor_diagonal(initial)

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