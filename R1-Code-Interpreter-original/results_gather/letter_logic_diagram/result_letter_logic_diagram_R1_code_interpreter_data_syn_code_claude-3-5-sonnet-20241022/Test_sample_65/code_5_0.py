def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def verify_initial_state(grid):
    # Check each row and column for duplicates
    for i in range(7):
        row_letters = set()
        col_letters = set()
        for j in range(7):
            if grid[i][j] != '':
                if grid[i][j] in row_letters:
                    return False
                row_letters.add(grid[i][j])
            if grid[j][i] != '':
                if grid[j][i] in col_letters:
                    return False
                col_letters.add(grid[j][i])
    return True

def is_valid(grid, row, col, letter, initial_grid):
    # Check initial constraints
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
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
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
            if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != letter:
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

    if not verify_initial_state(initial_grid):
        print("Invalid initial state")
        return False

    grid = [[''] * 7 for _ in range(7)]
    
    # Copy initial values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]

    def backtrack(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True

        # Skip pre-filled cells
        if initial_grid[row][col] != '':
            return backtrack(row, col + 1)

        # Try each letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    if backtrack():
        print_grid(grid)
        return True
    else:
        print("No solution found")
        return False

solve()