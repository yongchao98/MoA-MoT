def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['d', 'c', 'f', '', 'b', 'a', ''],
        ['', '', 'f', 'e', 'b', 'a', ''],
        ['c', '', '', 'e', 'b', 'a', ''],
        ['', '', 'b', 'a', 'g', 'd', ''],
        ['', '', '', 'g', 'd', 'c', ''],
        ['', '', '', 'd', '', '', 'e'],
        ['', '', '', 'c', 'f', '', '']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in minor_diagonal_indices if grid[i][j])
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter

    # Function to find missing letters in a list
    def find_missing_letters(lst):
        return list(all_letters - set(lst))

    # Fill the grid row by row
    for i in range(7):
        missing_letters = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that is not in the current column
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()