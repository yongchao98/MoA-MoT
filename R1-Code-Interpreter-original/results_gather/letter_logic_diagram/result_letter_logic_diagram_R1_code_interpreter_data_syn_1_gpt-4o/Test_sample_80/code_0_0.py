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
    # Check which letter can fit into all diagonal positions
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for letter in letters:
        if all(grid[r][c] in ('', letter) for r, c in diagonal_positions):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for r, c in diagonal_positions:
        grid[r][c] = diagonal_letter

    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[r])
                col_letters = set(grid[i][c] for i in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[r][c] = missing_letter

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(result)

solve_puzzle()