def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def forward_check(grid, row, col, letter):
        # Check if placing the letter leads to any immediate conflicts
        grid[row][col] = letter
        for i in range(7):
            if grid[row][i] == '' and not any(is_valid(grid, row, i, l) for l in 'abcdefg'):
                grid[row][col] = ''
                return False
            if grid[i][col] == '' and not any(is_valid(grid, i, col, l) for l in 'abcdefg'):
                grid[row][col] = ''
                return False
        grid[row][col] = ''
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter) and forward_check(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Determine the letter for the minor diagonal
    used_letters = set(grid[i][6-i] for i in range(7) if grid[i][6-i] != '')
    minor_diagonal_letter = next(letter for letter in 'abcdefg' if letter not in used_letters)

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    if solve(grid):
        return grid
    else:
        return None

# Initial grid with pre-filled values
grid = [
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

solution = solve_puzzle(grid)
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")