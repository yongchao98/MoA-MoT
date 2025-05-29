def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial):
    # FIRST: Must match initial configuration
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False

    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:  # Cell is on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def solve_puzzle():
    # Initial configuration
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]

    # Create working grid with initial values
    grid = [[''] * 7 for _ in range(7)]
    
    # First, copy all initial values and verify they don't violate constraints
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                if not is_valid(grid, i, j, initial[i][j], initial):
                    return None
                grid[i][j] = initial[i][j]

    def find_empty():
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
        # If on minor diagonal, try to match existing diagonal letter first
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    diagonal_letter = grid[i][j]
                    if is_valid(grid, row, col, diagonal_letter, initial):
                        grid[row][col] = diagonal_letter
                        if solve():
                            return True
                        grid[row][col] = ''
                    return False

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False

    if solve():
        print_solution(grid)
        return True
    return False

# Run the solver
solve_puzzle()