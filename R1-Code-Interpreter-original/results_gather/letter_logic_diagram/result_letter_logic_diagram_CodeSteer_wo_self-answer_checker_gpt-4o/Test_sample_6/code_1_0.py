def solve_puzzle():
    # Initial grid with empty cells represented as ''
    grid = [
        ['c', 'b', 'a', 'd', '', 'f', ''],
        ['b', '', '', 'e', '', '', ''],
        ['', 'd', '', 'f', 'g', '', ''],
        ['d', '', '', 'g', '', '', ''],
        ['', '', 'g', '', '', '', ''],
        ['', 'g', 'c', '', 'a', '', ''],
        ['', '', 'b', '', 'd', '', 'f']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is not present in any of the diagonal positions
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    used_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j]}
    diagonal_letter = (letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        missing_letters = list(letters - row_letters)
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that is not in the current column
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(f"<<<\n{result}\n>>>")

solve_puzzle()