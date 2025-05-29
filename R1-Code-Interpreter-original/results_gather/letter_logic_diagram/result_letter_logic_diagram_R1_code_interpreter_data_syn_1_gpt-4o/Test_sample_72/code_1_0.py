def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['e', 'c', 'b', '', 'f', 'd', 'g'],
        ['c', '', '', 'f', '', '', 'e'],
        ['b', 'a', '', '', 'g', 'e', 'c'],
        ['', 'f', '', '', 'e', 'c', ''],
        ['', 'd', '', 'e', '', '', ''],
        ['', '', 'e', 'c', 'b', '', 'f'],
        ['', '', 'c', '', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.discard(grid[i][6-i])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Track used letters in rows and columns
    used_in_row = [set(row) for row in grid]
    used_in_col = [set(grid[r][c] for r in range(7)) for c in range(7)]

    # Function to check if placing a letter is valid
    def is_valid(row, col, letter):
        return letter not in used_in_row[row] and letter not in used_in_col[col]

    # Backtracking function to fill the grid
    def fill_grid(row, col):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(row, col + 1)

        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                used_in_row[row].add(letter)
                used_in_col[col].add(letter)
                if fill_grid(row, col + 1):
                    return True
                grid[row][col] = ''
                used_in_row[row].remove(letter)
                used_in_col[col].remove(letter)

        return False

    # Start filling the grid
    fill_grid(0, 0)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()