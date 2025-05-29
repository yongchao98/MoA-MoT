def is_valid_move(grid, row, col, letter, initial):
    # Check if this position was pre-filled in initial grid
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
        # Find the required diagonal letter from initial grid
        diag_letter = None
        for i in range(7):
            if initial[i][6-i] != '':
                diag_letter = initial[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False

    return True

def verify_solution(grid, initial):
    # Verify initial constraints
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False

    # Verify rows and columns
    for i in range(7):
        row_set = set()
        col_set = set()
        for j in range(7):
            row_set.add(grid[i][j])
            col_set.add(grid[j][i])
        if len(row_set) != 7 or len(col_set) != 7:
            return False

    # Verify minor diagonal
    diag_set = set()
    for i in range(7):
        diag_set.add(grid[i][6-i])
    if len(diag_set) != 1:
        return False

    return True

def solve(grid, initial, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return verify_solution(grid, initial)

    # If this position is pre-filled, move to next position
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, row, col + 1)

    # If on minor diagonal, must use the correct letter
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if initial[i][6-i] != '':
                diag_letter = initial[i][6-i]
                break
        if diag_letter:
            if is_valid_move(grid, row, col, diag_letter, initial):
                grid[row][col] = diag_letter
                if solve(grid, initial, row, col + 1):
                    return True
            grid[row][col] = ''
            return False

    # Try each letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, row, col + 1):
                return True
            grid[row][col] = ''

    return False

# Initialize grids
initial = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

grid = [[''] * 7 for _ in range(7)]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')