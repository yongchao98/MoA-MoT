def solve_puzzle():
    # Initial grid with given letters
    grid = [
        ['', '', 'a', 'g', '', 'b', 'c'],
        ['', 'a', '', '', '', 'c', 'd'],
        ['a', 'g', 'e', 'b', '', '', ''],
        ['g', '', 'b', '', 'd', 'f', 'a'],
        ['e', 'b', '', '', 'f', '', 'g'],
        ['', '', 'd', 'f', 'a', '', ''],
        ['', 'd', 'f', 'a', 'g', 'e', 'b']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j])
    diagonal_letter = (letters - used_letters).pop()

    # Fill the minor diagonal with the determined letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = diagonal_letter

    # Track used letters in rows and columns
    row_used = [set(row) for row in grid]
    col_used = [set(grid[i][j] for i in range(7)) for j in range(7)]

    # Backtracking function to fill the grid
    def fill_grid(row, col):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(row, col + 1)

        for letter in letters:
            if letter not in row_used[row] and letter not in col_used[col]:
                grid[row][col] = letter
                row_used[row].add(letter)
                col_used[col].add(letter)

                if fill_grid(row, col + 1):
                    return True

                # Backtrack
                grid[row][col] = ''
                row_used[row].remove(letter)
                col_used[col].remove(letter)

        return False

    # Start filling the grid
    fill_grid(0, 0)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()