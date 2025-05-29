def solve_puzzle():
    # Initial grid with empty spaces represented by ''
    grid = [
        ['f', '', '', 'c', 'a', '', 'g'],
        ['e', 'b', 'c', 'a', 'd', 'g', ''],
        ['', 'c', '', '', 'g', 'f', 'e'],
        ['c', 'a', '', '', '', 'e', 'b'],
        ['a', 'd', 'g', '', 'e', '', 'c'],
        ['d', 'g', '', 'e', '', 'c', ''],
        ['', '', '', 'b', '', 'a', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][j] != letter for i, j in diagonal_positions):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Sets to track available letters for each row and column
    available_rows = [set('abcdefg') - set(row) for row in grid]
    available_cols = [set('abcdefg') - {grid[i][j] for i in range(7)} for j in range(7)]

    # Helper function to check if placing a letter is valid
    def is_valid(row, col, letter):
        # Check minor diagonal constraint
        if (row, col) in diagonal_positions and letter != diagonal_letter:
            return False
        return True

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in available_rows[row] & available_cols[col]:
            if is_valid(row, col, letter):
                grid[row][col] = letter
                available_rows[row].remove(letter)
                available_cols[col].remove(letter)

                if backtrack(row, col + 1):
                    return True

                # Backtrack
                grid[row][col] = ''
                available_rows[row].add(letter)
                available_cols[col].add(letter)

        return False

    # Start backtracking from the first cell
    if backtrack(0, 0):
        # Print the completed grid
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()