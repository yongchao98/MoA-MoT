def solve_puzzle():
    # Initial grid with empty cells represented as ''
    grid = [
        ['c', '', '', '', '', '', ''],
        ['', '', 'e', '', '', '', ''],
        ['d', 'e', 'g', '', 'f', 'c', 'b'],
        ['e', 'g', '', 'f', 'c', '', 'd'],
        ['', 'a', 'f', '', '', '', ''],
        ['', 'f', 'c', 'b', '', 'e', 'g'],
        ['f', 'c', '', '', 'e', 'g', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j]}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the determined letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Sets to track available letters for each row and column
    row_available = [set('abcdefg') - set(row) for row in grid]
    col_available = [set('abcdefg') - {grid[i][j] for i in range(7)} for j in range(7)]

    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in row_available[row] & col_available[col]:
            grid[row][col] = letter
            row_available[row].remove(letter)
            col_available[col].remove(letter)

            if backtrack(row, col + 1):
                return True

            # Backtrack
            grid[row][col] = ''
            row_available[row].add(letter)
            col_available[col].add(letter)

        return False

    if backtrack(0, 0):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()