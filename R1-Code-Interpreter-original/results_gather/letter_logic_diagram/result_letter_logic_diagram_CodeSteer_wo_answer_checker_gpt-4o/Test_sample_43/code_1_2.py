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

    # Initialize available letters for each row and column
    row_available = [set(all_letters) - set(row) for row in grid]
    col_available = [set(all_letters) - set(grid[i][j] for i in range(7)) for j in range(7)]

    def is_valid(row, col, letter):
        # Check if the letter is available in the row and column
        return letter in row_available[row] and letter in col_available[col]

    def solve():
        for row in range(7):
            for col in range(7):
                if grid[row][col] is None:
                    for letter in all_letters:
                        if is_valid(row, col, letter):
                            # Place the letter and update available sets
                            grid[row][col] = letter
                            row_available[row].remove(letter)
                            col_available[col].remove(letter)

                            if solve():
                                return True

                            # Backtrack
                            grid[row][col] = None
                            row_available[row].add(letter)
                            col_available[col].add(letter)
                    return False
        return True

    solve()

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()