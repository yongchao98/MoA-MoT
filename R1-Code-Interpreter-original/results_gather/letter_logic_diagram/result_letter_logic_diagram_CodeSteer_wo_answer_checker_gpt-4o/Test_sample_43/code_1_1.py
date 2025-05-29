def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['d', 'a', None, None, 'g', None, None],
        [None, None, None, 'g', 'e', None, None],
        ['f', None, 'g', None, 'c', None, None],
        [None, 'g', None, None, None, 'a', None],
        ['g', 'e', None, None, None, None, 'f'],
        ['e', None, None, None, None, 'b', None],
        [None, 'd', None, 'f', 'b', None, 'e']
    ]

    # Determine the letter for the minor diagonal
    diagonal_indices = [(i, 6-i) for i in range(7)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_indices if grid[i][j] is not None)
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check row and column
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid):
        for row in range(7):
            for col in range(7):
                if grid[row][col] is None:
                    for letter in all_letters:
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid):
                                return True
                            grid[row][col] = None
                    return False
        return True

    solve(grid)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()