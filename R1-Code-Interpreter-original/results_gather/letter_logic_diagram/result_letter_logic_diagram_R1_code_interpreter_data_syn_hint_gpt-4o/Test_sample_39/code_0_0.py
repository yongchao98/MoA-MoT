def fill_grid():
    # Initial grid setup
    grid = [
        ["", "b", "", "", "", "", ""],
        ["b", "a", "g", "d", "", "", "f"],
        ["", "", "d", "", "", "", "b"],
        ["g", "", "", "e", "", "", ""],
        ["", "", "", "", "b", "", ""],
        ["c", "", "f", "", "a", "", "d"],
        ["e", "f", "", "a", "g", "", ""]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    possible_letters = set("abcdefg")
    for i in range(7):
        for j in range(7):
            if grid[i][j] != "":
                possible_letters.discard(grid[i][j])

    # Choose a letter for the diagonal
    diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set("abcdefg")
        col_letters = set("abcdefg")
        for j in range(7):
            if grid[i][j] != "":
                row_letters.discard(grid[i][j])
            if grid[j][i] != "":
                col_letters.discard(grid[j][i])

        for j in range(7):
            if grid[i][j] == "":
                grid[i][j] = row_letters.pop()
            if grid[j][i] == "":
                grid[j][i] = col_letters.pop()

    # Print the filled grid
    for row in grid:
        print(",".join(row))

fill_grid()