def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, minor_diagonal_positions, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

    if col == 7:
        return solve(grid, minor_diagonal_letter, minor_diagonal_positions, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, minor_diagonal_positions, row, col + 1)

    for letter in 'abcdefg':
        if (row, col) in minor_diagonal_positions and letter != minor_diagonal_letter:
            continue

        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, minor_diagonal_positions, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_grid():
    grid = [
        ['a', 'f', 'b', 'd', '', 'c', ''],
        ['', '', 'd', '', 'c', '', 'a'],
        ['', '', 'e', '', '', '', 'f'],
        ['', '', 'c', '', '', 'f', ''],
        ['', '', '', 'a', '', 'b', ''],
        ['', '', '', 'f', 'b', '', ''],
        ['', 'a', 'f', '', '', 'e', '']
    ]

    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    for minor_diagonal_letter in 'abcdefg':
        if solve(grid, minor_diagonal_letter, minor_diagonal_positions):
            for row in grid:
                print(','.join(row))
            return

fill_grid()