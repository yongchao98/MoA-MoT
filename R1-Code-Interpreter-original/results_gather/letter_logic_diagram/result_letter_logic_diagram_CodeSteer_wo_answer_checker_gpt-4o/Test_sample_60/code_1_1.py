def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        [None, 'c', 'g', 'd', 'b', None, None],
        [None, None, 'd', 'b', None, None, None],
        [None, 'd', None, None, None, None, None],
        [None, 'b', None, None, None, 'c', None],
        ['b', 'e', None, 'a', 'c', None, 'd'],
        ['e', None, None, None, None, 'd', None],
        ['f', None, None, 'g', 'd', None, None]
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None)
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] is None:
                    for letter in all_letters:
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            if solve(grid):
                                return True
                            grid[i][j] = None
                    return False
        return True

    solve(grid)

    # Print the result
    for row in grid:
        print(','.join(row))

solve_puzzle()