def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['d', 'f', 'e', 'b', 'a', 'g', ''],
        ['f', '', 'b', 'a', 'g', '', ''],
        ['e', 'b', '', '', 'c', 'd', ''],
        ['b', '', 'g', 'c', '', 'f', ''],
        ['b', 'a', '', 'c', 'd', '', ''],
        ['a', '', '', 'd', 'f', 'e', 'b'],
        ['', '', '', '', 'e', '', 'a']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            letters.discard(grid[r][c])
    minor_diagonal_letter = letters.pop()

    # Fill the minor diagonal with the chosen letter
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter

    # Fill the rest of the grid
    for r in range(7):
        row_letters = set('abcdefg') - set(grid[r])
        for c in range(7):
            if grid[r][c] == '':
                for letter in row_letters:
                    # Check if the letter can be placed in this column
                    if all(grid[i][c] != letter for i in range(7)):
                        grid[r][c] = letter
                        row_letters.remove(letter)
                        break

    # Print the final grid
    for row in grid:
        print(','.join(row))

solve_puzzle()