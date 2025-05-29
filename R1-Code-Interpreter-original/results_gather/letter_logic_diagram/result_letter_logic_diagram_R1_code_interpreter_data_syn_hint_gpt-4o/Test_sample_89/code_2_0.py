def fill_grid():
    # Initial grid with given letters
    grid = [
        ['', '', 'a', 'g', '', 'b', 'c'],
        ['', 'a', '', '', '', 'c', 'd'],
        ['a', 'g', 'e', 'b', '', '', ''],
        ['g', '', 'b', '', 'd', 'f', 'a'],
        ['e', 'b', '', '', 'f', '', 'g'],
        ['', '', 'd', 'f', 'a', '', ''],
        ['', 'd', 'f', 'a', 'g', 'e', 'b']
    ]

    # Function to check if a letter can be placed on the minor diagonal
    def can_place_on_diagonal(letter):
        for i in range(7):
            if letter in grid[i] or letter in [grid[k][6-i] for k in range(7)]:
                return False
        return True

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if can_place_on_diagonal(letter):
            minor_diagonal_letter = letter
            break

    # If no suitable letter is found, raise an error
    if minor_diagonal_letter is None:
        raise ValueError("No suitable letter found for the minor diagonal.")

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can fit in this position
                col_letters = set(grid[k][j] for k in range(7))
                possible_letter = (row_letters - col_letters).pop()
                grid[i][j] = possible_letter
                row_letters.remove(possible_letter)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

fill_grid()