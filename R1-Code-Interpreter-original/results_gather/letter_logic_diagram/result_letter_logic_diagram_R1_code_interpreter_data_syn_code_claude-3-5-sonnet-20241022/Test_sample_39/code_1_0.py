def verify_fixed_cells(grid, fixed):
    for i in range(7):
        for j in range(7):
            if fixed[i][j] != '' and grid[i][j] != fixed[i][j]:
                return False
    return True

def is_valid_move(grid, row, col, letter, fixed):
    # Must match fixed cell if it exists
    if fixed[row][col] != '' and fixed[row][col] != letter:
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

def solve_grid(fixed):
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # First, place all fixed values
    for i in range(7):
        for j in range(7):
            if fixed[i][j] != '':
                grid[i][j] = fixed[i][j]

    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return verify_fixed_cells(grid, fixed)

        row = pos // 7
        col = pos % 7

        # Skip if cell is fixed
        if fixed[row][col] != '':
            return backtrack(pos + 1)

        # Try each letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, fixed):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''

        return False

    if backtrack():
        return grid
    return None

# Initial configuration
fixed = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

print('<<<')
result = solve_grid(fixed)
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')