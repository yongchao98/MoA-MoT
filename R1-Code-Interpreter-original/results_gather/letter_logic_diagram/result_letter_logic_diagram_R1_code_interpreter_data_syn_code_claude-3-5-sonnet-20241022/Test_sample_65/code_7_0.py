def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def validate_initial_grid(initial_grid):
    # Check for any obvious conflicts in initial values
    for i in range(7):
        row_letters = set()
        col_letters = set()
        for j in range(7):
            if initial_grid[i][j] != '':
                if initial_grid[i][j] in row_letters:
                    return False
                row_letters.add(initial_grid[i][j])
            if initial_grid[j][i] != '':
                if initial_grid[j][i] in col_letters:
                    return False
                col_letters.add(initial_grid[j][i])
    return True

def is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
    # Must match initial value if present
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Must match minor diagonal letter if on minor diagonal
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

    if not validate_initial_grid(initial_grid):
        print("Invalid initial state")
        return False

    # Find minor diagonal letter from initial values
    minor_diag_letter = 'f'  # We can see 'f' appears on the minor diagonal in initial grid

    grid = [[''] * 7 for _ in range(7)]
    
    # Fill in initial values and minor diagonal
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]
            elif i + j == 6:  # Minor diagonal
                grid[i][j] = minor_diag_letter

    def backtrack(pos=0):
        if pos == 49:
            return True

        row = pos // 7
        col = pos % 7

        # Skip pre-filled cells
        if grid[row][col] != '':
            return backtrack(pos + 1)

        # Try each possible letter
        letters = list('abcdefg')
        # Prioritize minor_diag_letter if on minor diagonal
        if row + col == 6:
            letters.remove(minor_diag_letter)
            letters.insert(0, minor_diag_letter)

        for letter in letters:
            if is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
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