def solve_puzzle():
    # Initial grid with given letters and empty spaces
    grid = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diagonal_letter = grid[i][6-i]
            break

    # If no letter is found on the diagonal, choose one (e.g., 'a')
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Sets to track used letters in each row and column
    row_sets = [set(row) for row in grid]
    col_sets = [set(grid[i][j] for i in range(7)) for j in range(7)]

    # Function to solve the grid using backtracking
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if letter not in row_sets[row] and letter not in col_sets[col]:
                grid[row][col] = letter
                row_sets[row].add(letter)
                col_sets[col].add(letter)

                if backtrack(row, col + 1):
                    return True

                # Backtrack
                grid[row][col] = ''
                row_sets[row].remove(letter)
                col_sets[col].remove(letter)

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()