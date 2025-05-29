def solve_puzzle():
    # Initial grid setup
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
    # Check which letter can fit into all diagonal positions
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letters.discard(grid[i][j])

    # Choose a letter for the diagonal
    diagonal_letter = None
    for letter in letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                col_letters = set(grid[k][j] for k in range(7) if grid[k][j] != '')
                possible_letters = row_letters - col_letters
                if possible_letters:
                    grid[i][j] = possible_letters.pop()
                else:
                    # If no possible letter is found, print an error message
                    print(f"Error: No valid letter for cell ({i}, {j})")
                    return

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()