def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, num, original):
    # Must match original if cell was pre-filled
    if original[row][col] != '' and num != original[row][col]:
        return False

    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False

    # Check column    
    for x in range(7):
        if grid[x][col] == num:
            return False

    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != num:
                return False

    return True

# Original puzzle
original = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

# Initialize solution grid
grid = [['' for _ in range(7)] for _ in range(7)]

# First, copy all pre-filled values
for i in range(7):
    for j in range(7):
        if original[i][j] != '':
            grid[i][j] = original[i][j]

def solve_grid():
    # Find minor diagonal letter (must be consistent)
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if not diag_letter:
        # If no diagonal letter found, try each possibility
        for letter in 'abcdefg':
            valid = True
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    valid = False
                    break
            if valid:
                diag_letter = letter
                break

    def find_empty():
        # First fill minor diagonal
        for i in range(7):
            if i + (6-i) == 6 and grid[i][6-i] == '':
                return (i, 6-i)
        # Then fill other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None

    def solve():
        pos = find_empty()
        if not pos:
            return True

        row, col = pos

        # If on minor diagonal, only try diagonal letter
        if row + col == 6:
            if is_valid_move(grid, row, col, diag_letter, original):
                grid[row][col] = diag_letter
                if solve():
                    return True
                grid[row][col] = ''
            return False

        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, original):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''

        return False

    if solve():
        print_grid(grid)
        return True
    return False

solve_grid()