def verify_initial_values(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def verify_rows_cols(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set()
        col_set = set()
        for j in range(7):
            if grid[i][j] in row_set or grid[j][i] in col_set:
                return False
            row_set.add(grid[i][j])
            col_set.add(grid[j][i])
    return True

def verify_minor_diagonal(grid):
    # Check minor diagonal (top-right to bottom-left)
    diagonal_letter = None
    for i in range(7):
        j = 6 - i
        if diagonal_letter is None:
            diagonal_letter = grid[i][j]
        elif grid[i][j] != diagonal_letter:
            return False
    return True

def is_valid_move(grid, row, col, letter, initial_grid):
    # Check initial constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Temporarily place the letter
    temp = grid[row][col]
    grid[row][col] = letter

    # Check row
    row_letters = set()
    for j in range(7):
        if grid[row][j] != '':
            if grid[row][j] in row_letters:
                grid[row][col] = temp
                return False
            row_letters.add(grid[row][j])

    # Check column
    col_letters = set()
    for i in range(7):
        if grid[i][col] != '':
            if grid[i][col] in col_letters:
                grid[row][col] = temp
                return False
            col_letters.add(grid[i][col])

    # Check minor diagonal
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                if diagonal_letter is None:
                    diagonal_letter = grid[i][j]
                elif grid[i][j] != diagonal_letter:
                    grid[row][col] = temp
                    return False

    grid[row][col] = temp
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:
        return verify_rows_cols(grid) and verify_minor_diagonal(grid) and verify_initial_values(grid, initial_grid)

    row = pos // 7
    col = pos % 7

    # If cell is pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)

    # For minor diagonal cells, try to use existing diagonal letter first
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter and is_valid_move(grid, row, col, diagonal_letter, initial_grid):
            grid[row][col] = diagonal_letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''

    # Try each letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''

    return False

initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

grid = [[''] * 7 for _ in range(7)]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")