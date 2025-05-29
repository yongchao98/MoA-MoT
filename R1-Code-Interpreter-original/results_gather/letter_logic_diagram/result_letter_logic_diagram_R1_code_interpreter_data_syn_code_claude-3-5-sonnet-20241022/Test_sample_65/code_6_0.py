def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def get_minor_diagonal_value(initial_grid):
    # First, find any pre-filled value on the minor diagonal
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
    # Check initial value constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check minor diagonal constraint
    if row + col == 6 and letter != minor_diag_letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    return True

def solve():
    initial_grid = [
        ['a','','d','','b','g','c'],
        ['f','d','','','','','a'],
        ['d','','','','','',''],
        ['','','','','','f',''],
        ['','g','','a','f','','e'],
        ['','c','','f','','e',''],
        ['','a','','','e','','']
    ]

    # Get the minor diagonal letter
    minor_diag_letter = get_minor_diagonal_value(initial_grid)
    if not minor_diag_letter:
        # If no minor diagonal letter is pre-filled, we need to find one that works
        for test_letter in 'abcdefg':
            if try_solve_with_minor_diag(initial_grid, test_letter):
                return True
        return False
    else:
        return try_solve_with_minor_diag(initial_grid, minor_diag_letter)

def try_solve_with_minor_diag(initial_grid, minor_diag_letter):
    grid = [[''] * 7 for _ in range(7)]
    
    # First fill in all initial values and minor diagonal
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]
            elif i + j == 6:  # Minor diagonal
                grid[i][j] = minor_diag_letter

    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True

        row = pos // 7
        col = pos % 7

        # Skip pre-filled cells and minor diagonal
        if grid[row][col] != '':
            return backtrack(pos + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''

        return False

    # Verify initial state is valid before proceeding
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                temp = grid[i][j]
                grid[i][j] = ''
                if not is_valid(grid, i, j, temp, initial_grid, minor_diag_letter):
                    return False
                grid[i][j] = temp

    if backtrack():
        print_grid(grid)
        return True
    return False

solve()