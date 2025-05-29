def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

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
    # Initialize the grid with given values
    initial_grid = [
        ['a','','d','','b','g','c'],
        ['f','d','','','','','a'],
        ['d','','','','','',''],
        ['','','','','','f',''],
        ['','g','','a','f','','e'],
        ['','c','','f','','e',''],
        ['','a','','','e','','']
    ]

    # We can see that 'f' must be on the minor diagonal because of position [4,5]
    minor_diag_letter = 'f'

    def try_solve():
        grid = [[''] * 7 for _ in range(7)]
        
        # First, fill in all initial values
        for i in range(7):
            for j in range(7):
                if initial_grid[i][j] != '':
                    grid[i][j] = initial_grid[i][j]
                elif i + j == 6:  # Fill minor diagonal with 'f'
                    grid[i][j] = minor_diag_letter

        def backtrack(row=0, col=0):
            if col == 7:
                row += 1
                col = 0
            if row == 7:
                return True

            # If cell is already filled, move to next cell
            if grid[row][col] != '':
                return backtrack(row, col + 1)

            # Try each possible letter
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
                    grid[row][col] = letter
                    if backtrack(row, col + 1):
                        return True
                    grid[row][col] = ''

            return False

        # First verify that initial state is valid
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

        if backtrack(0, 0):
            print_grid(grid)
            return True
        return False

    if not try_solve():
        print("No solution found")

solve()