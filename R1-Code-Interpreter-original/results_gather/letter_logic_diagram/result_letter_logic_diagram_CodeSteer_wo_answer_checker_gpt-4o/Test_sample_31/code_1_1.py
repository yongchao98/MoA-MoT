def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['e', 'a', None, None, None, 'c', 'f'],
        [None, 'b', None, None, 'c', 'f', None],
        ['b', 'g', 'd', 'c', 'f', None, 'a'],
        ['g', 'd', None, 'f', 'e', None, 'b'],
        ['d', None, None, 'e', None, None, 'g'],
        [None, 'f', 'e', None, 'b', None, 'd'],
        [None, None, 'a', None, 'g', 'd', 'c']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = 'e'
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            # If we've filled all rows, print the solution
            for r in grid:
                print(','.join(r))
            return True

        if col == 7:
            # Move to the next row
            return solve(grid, row + 1, 0)

        if grid[row][col] is not None:
            # Skip pre-filled cells
            return solve(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = None  # Backtrack

        return False

    solve(grid)

solve_puzzle()