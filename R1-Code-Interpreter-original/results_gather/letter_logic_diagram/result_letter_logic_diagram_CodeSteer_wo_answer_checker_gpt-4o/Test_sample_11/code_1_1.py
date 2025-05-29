def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        [None, None, 'd', 'g', None, 'b', 'f'],
        [None, 'd', None, None, None, 'f', None],
        ['d', 'g', None, None, 'f', None, None],
        ['g', None, 'b', 'f', 'e', 'c', 'd'],
        ['a', 'b', 'f', None, 'c', 'd', None],
        [None, None, None, 'c', 'd', 'g', None],
        ['f', None, 'c', 'd', 'g', 'a', 'b']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Helper function to check if placing a letter is valid
    def is_valid(i, j, letter):
        # Check row and column
        if letter in grid[i] or letter in (grid[x][j] for x in range(7)):
            return False
        return True

    # Backtracking function
    def backtrack(cell=0):
        if cell == 49:  # All cells are filled
            return True
        i, j = divmod(cell, 7)
        if grid[i][j] is not None:  # Skip pre-filled cells
            return backtrack(cell + 1)
        for letter in all_letters:
            if is_valid(i, j, letter):
                grid[i][j] = letter
                if backtrack(cell + 1):
                    return True
                grid[i][j] = None  # Backtrack
        return False

    # Start backtracking from the first cell
    backtrack()

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()