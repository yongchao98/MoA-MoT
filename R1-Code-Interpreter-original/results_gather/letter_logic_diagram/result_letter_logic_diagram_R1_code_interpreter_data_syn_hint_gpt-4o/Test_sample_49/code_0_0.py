def solve_puzzle():
    # Initial grid with given letters
    grid = [
        ["", "", "", "", "d", "", "e"],
        ["c", "f", "", "", "", "g", ""],
        ["f", "b", "", "", "", "a", "c"],
        ["", "d", "", "e", "", "c", ""],
        ["", "g", "", "a", "", "", "b"],
        ["", "", "a", "c", "", "", ""],
        ["", "", "", "", "b", "d", ""]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    letters = set("abcdefg")
    for letter in letters:
        if all(grid[i][6-i] == "" or grid[i][6-i] == letter for i in range(7)):
            minor_diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        missing_letters = letters - row_letters
        for j in range(7):
            if grid[i][j] == "":
                # Find a letter that can be placed in this position
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):  # Check column constraint
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()