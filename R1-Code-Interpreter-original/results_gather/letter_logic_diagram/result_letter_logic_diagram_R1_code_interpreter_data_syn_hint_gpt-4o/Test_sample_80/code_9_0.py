def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue
        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''

    return False

def fill_grid():
    grid = [
        ['', '', '', '', '', 'a', ''],
        ['', 'e', '', 'd', '', '', ''],
        ['', 'c', '', 'a', '', 'g', 'f'],
        ['', '', '', 'b', '', 'f', 'e'],
        ['d', 'a', '', '', '', 'e', ''],
        ['', '', '', 'f', 'e', '', ''],
        ['b', 'g', 'f', 'e', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[r][c] in ('', letter) for r, c in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]):
            minor_diagonal_letter = letter
            break

    if minor_diagonal_letter is None:
        return "No valid solution for the minor diagonal"

    # Solve the grid
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

fill_grid()