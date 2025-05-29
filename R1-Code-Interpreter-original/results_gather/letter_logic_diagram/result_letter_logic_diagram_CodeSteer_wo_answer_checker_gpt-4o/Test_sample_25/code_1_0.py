def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['f', None, 'd', None, None, 'g', 'a'],
        [None, None, 'e', None, None, None, None],
        [None, 'e', None, 'g', None, None, 'c'],
        [None, None, None, 'a', None, 'c', 'd'],
        [None, 'g', None, 'f', 'c', 'd', 'e'],
        [None, None, None, None, None, None, None],
        [None, None, 'c', 'd', 'e', None, None]
    ]

    all_letters = set('abcdefg')
    n = 7

    # Determine possible letters for the minor diagonal
    diagonal_positions = [(i, n - 1 - i) for i in range(n)]
    used_diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None)
    possible_diagonal_letters = all_letters - used_diagonal_letters

    # Helper function to check if placing a letter is valid
    def is_valid(row, col, letter):
        # Check row and column
        for i in range(n):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    # Backtracking function with logging
    def backtrack(row, col):
        if row == n:
            return True
        if col == n:
            return backtrack(row + 1, 0)
        if grid[row][col] is not None:
            return backtrack(row, col + 1)

        for letter in all_letters:
            if is_valid(row, col, letter):
                grid[row][col] = letter
                print(f"Placing {letter} at ({row}, {col})")
                if backtrack(row, col + 1):
                    return True
                print(f"Backtracking from ({row}, {col})")
                grid[row][col] = None

        return False

    # Try each possible letter for the minor diagonal
    for diagonal_letter in possible_diagonal_letters:
        # Fill the minor diagonal with the chosen letter
        for i, j in diagonal_positions:
            grid[i][j] = diagonal_letter

        # Start backtracking from the first cell
        if backtrack(0, 0):
            # Verify the solution
            def verify_solution():
                for i in range(n):
                    if set(grid[i]) != all_letters:
                        return False
                    if set(grid[j][i] for j in range(n)) != all_letters:
                        return False
                if all(grid[i][n - 1 - i] == diagonal_letter for i in range(n)):
                    return True
                return False

            if verify_solution():
                # Print the completed grid
                for row in grid:
                    print(','.join(row))
                return

    print("No valid solution found.")

solve_puzzle()