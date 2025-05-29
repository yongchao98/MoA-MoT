def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row, col):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6 - i] in ('', letter) for i in range(7)):
            minor_diagonal_letter = letter
            break

    if minor_diagonal_letter is None:
        print("No valid letter found for the minor diagonal.")
        return

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    if solve(grid, 0, 0):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid setup
grid = [
    ['', 'a', 'e', 'g', 'f', 'b', ''],
    ['a', '', 'g', 'f', '', '', ''],
    ['', '', '', '', 'c', 'd', 'a'],
    ['', 'f', 'b', 'c', '', 'a', ''],
    ['', 'b', 'c', '', 'a', '', 'g'],
    ['b', '', '', '', '', 'g', ''],
    ['c', '', '', 'e', 'g', '', '']
]

solve_puzzle(grid)