def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, original):
    # Check if this position was pre-filled in original
    if original[row][col] != '' and letter != original[row][col]:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False

    # Check minor diagonal
    if row + col == 6:  # If on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def solve():
    # Original puzzle with pre-filled values
    original = [
        ['a','e','','f','g','c',''],
        ['','b','','','','d','a'],
        ['b','f','g','c','d','a','e'],
        ['f','','c','','a','','b'],
        ['g','c','d','a','e','b','f'],
        ['c','','','e','b','','g'],
        ['d','','','','f','g','c']
    ]

    # Working grid
    grid = [row[:] for row in original]

    def find_next_empty():
        # First fill minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j)
        # Then fill other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None

    def solve_recursive():
        cell = find_next_empty()
        if not cell:
            return True

        row, col = cell
        possible_letters = 'abcdefg'

        # If on minor diagonal, use existing diagonal letter if any exists
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    possible_letters = grid[i][j]
                    break

        for letter in possible_letters:
            if is_valid(grid, row, col, letter, original):
                grid[row][col] = letter
                if solve_recursive():
                    return True
                grid[row][col] = ''

        return False

    if solve_recursive():
        print_grid(grid)
        return True
    else:
        print("No solution exists")
        return False

solve()