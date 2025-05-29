def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != diagonal_letter) or not is_valid(grid, row, col, letter):
            continue
        grid[row][col] = letter
        print(f"Placing {letter} at ({row}, {col})")  # Debugging output
        if solve(grid, diagonal_letter, row, col + 1):
            return True
        print(f"Backtracking from ({row}, {col})")  # Debugging output
        grid[row][col] = ''

    return False

def solve_puzzle():
    grid = [
        ['', '', 'e', '', '', '', ''],
        ['g', 'e', 'f', '', '', '', 'c'],
        ['', '', '', 'a', '', '', ''],
        ['', '', '', '', '', '', 'e'],
        ['', 'a', 'd', 'c', '', '', ''],
        ['', 'd', 'c', 'g', '', '', 'b'],
        ['d', 'c', 'g', '', 'f', 'b', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_indices if grid[i][j]}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Solve the puzzle using backtracking
    if solve(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()