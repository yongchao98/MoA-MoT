def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        [None, None, None, None, None, 'a', None],
        [None, 'e', None, 'd', None, None, None],
        [None, 'c', None, 'a', None, 'g', 'f'],
        [None, None, None, 'b', None, 'f', 'e'],
        ['d', 'a', None, None, None, 'e', None],
        [None, None, None, 'f', 'e', None, None],
        ['b', 'g', 'f', 'e', None, None, None]
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None)
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - used_letters).pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def backtrack(grid, row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(grid, row + 1, 0)
        if grid[row][col] is not None:
            return backtrack(grid, row, col + 1)

        for letter in all_letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(grid, row, col + 1):
                    return True
                grid[row][col] = None

        return False

    # Start backtracking from the first cell
    backtrack(grid, 0, 0)

    # Validate the grid
    for i, j in diagonal_positions:
        assert grid[i][j] == diagonal_letter, "Diagonal constraint violated"

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()