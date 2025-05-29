def solve_puzzle():
    # Initial grid setup
    grid = [
        ["", "e", "b", "", "", "f", ""],
        ["", "", "d", "c", "f", "g", "a"],
        ["", "", "c", "", "", "", "e"],
        ["", "c", "f", "g", "a", "e", ""],
        ["c", "", "g", "", "", "", "d"],
        ["f", "g", "a", "e", "b", "d", "c"],
        ["", "", "e", "b", "", "c", "f"]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is missing from the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set("abcdefg")
    minor_diagonal_letters = set(grid[i][j] for i, j in minor_diagonal_indices if grid[i][j])
    missing_letter = list(letters - minor_diagonal_letters)[0]

    # Fill the minor diagonal with the missing letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = missing_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        missing_row_letters = list(letters - row_letters)
        for j in range(7):
            if grid[i][j] == "":
                # Find a letter that is not in the current column
                for letter in missing_row_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_row_letters.remove(letter)
                        break

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()