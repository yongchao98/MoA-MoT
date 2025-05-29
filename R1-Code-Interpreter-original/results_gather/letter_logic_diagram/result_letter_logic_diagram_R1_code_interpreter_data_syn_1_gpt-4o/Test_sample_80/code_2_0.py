def solve_puzzle():
    # Initial grid setup
    grid = [
        ['', '', '', '', '', 'a', ''],
        ['', 'e', '', 'd', '', '', ''],
        ['', 'c', '', 'a', '', 'g', 'f'],
        ['', '', '', 'b', '', 'f', 'e'],
        ['d', 'a', '', '', '', 'e', ''],
        ['', '', '', 'f', 'e', '', ''],
        ['b', 'g', 'f', 'e', '', '', '']
    ]

    # Possible letters
    letters = set('abcdefg')

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for letter in letters:
        if all(grid[r][c] in ('', letter) for r, c in diagonal_positions):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for r, c in diagonal_positions:
        grid[r][c] = diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(r, c):
        if r == 7:
            return True
        if c == 7:
            return backtrack(r + 1, 0)
        if grid[r][c] != '':
            return backtrack(r, c + 1)

        row_letters = set(grid[r])
        col_letters = set(grid[i][c] for i in range(7))
        possible_letters = letters - row_letters - col_letters

        for letter in possible_letters:
            grid[r][c] = letter
            if backtrack(r, c + 1):
                return True
            grid[r][c] = ''

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(result)

solve_puzzle()