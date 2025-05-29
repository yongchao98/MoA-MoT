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
    def fill_grid():
        # Find the first empty cell
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    # Try placing each possible letter
                    for letter in possible_values[i][j]:
                        grid[i][j] = letter
                        # Save current state of possible values
                        original_row = [possible_values[i][k].copy() for k in range(7)]
                        original_col = [possible_values[k][j].copy() for k in range(7)]
                        # Update possible values
                        for k in range(7):
                            possible_values[i][k].discard(letter)
                            possible_values[k][j].discard(letter)
                        # Recursively fill the rest of the grid
                        if fill_grid():
                            return True
                        # Backtrack
                        grid[i][j] = ''
                        for k in range(7):
                            possible_values[i][k] = original_row[k]
                            possible_values[k][j] = original_col[k]
                    return False
        return True

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