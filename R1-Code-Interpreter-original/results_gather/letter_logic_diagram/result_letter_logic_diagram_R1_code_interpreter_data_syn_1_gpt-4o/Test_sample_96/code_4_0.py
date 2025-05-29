def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(n):
        if grid[i][n - 1 - i] != '':
            if minor_diagonal_letter is None:
                minor_diagonal_letter = grid[i][n - 1 - i]
            elif minor_diagonal_letter != grid[i][n - 1 - i]:
                raise ValueError("Conflicting letters on the minor diagonal")

    # If no letter is pre-filled on the minor diagonal, choose one
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Arbitrarily choose 'a'

    # Fill the minor diagonal
    for i in range(n):
        grid[i][n - 1 - i] = minor_diagonal_letter

    # Track used letters in rows and columns
    row_used = [set(row) for row in grid]
    col_used = [set(grid[i][j] for i in range(n)) for j in range(n)]

    def backtrack(r, c):
        if r == n:
            return True
        if c == n:
            return backtrack(r + 1, 0)
        if grid[r][c] != '':
            return backtrack(r, c + 1)

        for letter in letters:
            if letter not in row_used[r] and letter not in col_used[c]:
                grid[r][c] = letter
                row_used[r].add(letter)
                col_used[c].add(letter)

                if backtrack(r, c + 1):
                    return True

                # Backtrack
                grid[r][c] = ''
                row_used[r].remove(letter)
                col_used[c].remove(letter)

        return False

    if not backtrack(0, 0):
        raise ValueError("No solution found")

    return grid

# Initial grid
grid = [
    ['g', '', '', '', '', '', ''],
    ['', '', 'c', '', '', 'b', 'g'],
    ['d', '', '', '', '', '', 'g'],
    ['c', '', '', '', '', 'f', ''],
    ['e', '', 'b', 'g', 'f', '', 'c'],
    ['a', 'b', '', 'f', '', '', ''],
    ['', '', 'f', '', 'c', '', 'a']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))