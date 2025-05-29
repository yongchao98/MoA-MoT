def fill_grid():
    # Initial grid with given letters
    grid = [
        ["", "", "", "b", "", "e", ""],
        ["f", "", "", "d", "", "c", "a"],
        ["", "", "", "e", "", "a", "f"],
        ["b", "d", "", "", "", "f", ""],
        ["d", "e", "", "", "", "", ""],
        ["", "", "", "f", "g", "", ""],
        ["c", "", "f", "", "b", "", ""]
    ]

    # Determine the letter for the minor diagonal
    # We will try each letter from 'a' to 'g' and see which one fits
    letters = set("abcdefg")
    for letter in letters:
        # Check if the letter can be placed on the minor diagonal
        if all(grid[i][6-i] == "" or grid[i][6-i] == letter for i in range(7)):
            # Try to fill the grid with this letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter

            # Now fill the rest of the grid
            for i in range(7):
                for j in range(7):
                    if grid[i][j] == "":
                        # Find a letter that can be placed here
                        row_letters = set(grid[i])
                        col_letters = set(grid[k][j] for k in range(7))
                        possible_letters = letters - row_letters - col_letters
                        if possible_letters:
                            grid[i][j] = possible_letters.pop()

            # Check if the grid is valid
            if all(set(row) == letters for row in grid) and all(set(grid[i][j] for i in range(7)) == letters for j in range(7)):
                return grid

    return None

# Fill the grid and print the result
filled_grid = fill_grid()
if filled_grid:
    for row in filled_grid:
        print(",".join(row))