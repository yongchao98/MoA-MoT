def solve_puzzle(grid):
    # Initialize sets to track available letters for each row and column
    letters = set('abcdefg')
    row_avail = [letters - set(row) for row in grid]
    col_avail = [letters - set(grid[i][j] for i in range(7)) for j in range(7)]

    # Choose a letter for the minor diagonal
    diagonal_letter = 'e'  # This can be any letter that fits the constraints
    for i in range(7):
        grid[i][6 - i] = diagonal_letter
        if diagonal_letter in row_avail[i]:
            row_avail[i].remove(diagonal_letter)
        if diagonal_letter in col_avail[6 - i]:
            col_avail[6 - i].remove(diagonal_letter)

    def backtrack(row, col):
        if row == 7:
            return True  # All rows are filled
        if col == 7:
            return backtrack(row + 1, 0)  # Move to the next row

        if grid[row][col] != '':
            return backtrack(row, col + 1)  # Skip pre-filled cells

        for letter in row_avail[row] & col_avail[col]:
            grid[row][col] = letter
            row_avail[row].remove(letter)
            col_avail[col].remove(letter)

            if backtrack(row, col + 1):
                return True

            # Undo the move
            grid[row][col] = ''
            row_avail[row].add(letter)
            col_avail[col].add(letter)

        return False

    backtrack(0, 0)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid with empty cells as ''
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

print(solve_puzzle(grid))