def solve_puzzle(grid):
    from collections import defaultdict

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Pre-fill the grid and update possible values
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letter = grid[i][j]
                possible_values[i][j] = {letter}
                for k in range(7):
                    possible_values[i][k].discard(letter)
                    possible_values[k][j].discard(letter)

    # Choose a common letter for the minor diagonal
    minor_diagonal = [(i, 6 - i) for i in range(7)]
    common_letter = None
    for letter in 'abcdefg':
        if all(grid[i][j] == '' or grid[i][j] == letter for i, j in minor_diagonal):
            common_letter = letter
            break

    if common_letter is None:
        print("No valid common letter for the minor diagonal")
        return

    for i, j in minor_diagonal:
        grid[i][j] = common_letter
        possible_values[i][j] = {common_letter}
        for k in range(7):
            possible_values[i][k].discard(common_letter)
            possible_values[k][j].discard(common_letter)

    # Backtracking function to fill the grid
    def fill_grid(row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(row, col + 1)

        for letter in possible_values[row][col]:
            grid[row][col] = letter
            # Update possible values
            original_values = [possible_values[row][k] for k in range(7)]
            for k in range(7):
                possible_values[row][k].discard(letter)
                possible_values[k][col].discard(letter)

            if fill_grid(row, col + 1):
                return True

            # Backtrack
            grid[row][col] = ''
            for k in range(7):
                possible_values[row][k] = original_values[k]

        return False

    if fill_grid():
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid
grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['a', '', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

solve_puzzle(grid)