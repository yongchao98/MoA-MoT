def solve_puzzle():
    # Initial grid with empty spaces represented by ''
    grid = [
        ['g', '', '', 'c', '', '', 'd'],
        ['', 'b', '', 'e', '', 'd', ''],
        ['b', '', 'e', 'f', 'd', '', 'a'],
        ['c', 'e', 'f', '', '', 'a', 'b'],
        ['e', 'f', 'd', 'g', 'a', 'b', ''],
        ['f', '', '', 'a', 'b', '', ''],
        ['', 'g', 'a', 'b', '', '', 'f']
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.intersection_update(grid[i][6-i])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Sets to track available letters for each row and column
    row_sets = [set('abcdefg') - set(row) for row in grid]
    col_sets = [set('abcdefg') - set(grid[i][j] for i in range(7)) for j in range(7)]

    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in row_sets[row] & col_sets[col]:
            grid[row][col] = letter
            row_sets[row].remove(letter)
            col_sets[col].remove(letter)

            if backtrack(row, col + 1):
                return True

            # Backtrack
            grid[row][col] = ''
            row_sets[row].add(letter)
            col_sets[col].add(letter)

        return False

    backtrack(0, 0)

    # Format the output
    for row in grid:
        print(','.join(row))

solve_puzzle()