def fill_grid():
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

    # Determine the letter for the minor diagonal
    # We will try each letter from 'a' to 'g' and see if it fits
    letters = set('abcdefg')
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Check which letter can be placed on the minor diagonal
    for letter in letters:
        can_place = True
        for r, c in minor_diagonal_positions:
            if grid[r][c] != '' and grid[r][c] != letter:
                can_place = False
                break
            # Check row and column constraints
            if letter in grid[r] or any(grid[i][c] == letter for i in range(7)):
                can_place = False
                break
        if can_place:
            for r, c in minor_diagonal_positions:
                grid[r][c] = letter
            break

    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in letters:
                    if letter not in grid[r] and all(grid[i][c] != letter for i in range(7)):
                        grid[r][c] = letter
                        break

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()