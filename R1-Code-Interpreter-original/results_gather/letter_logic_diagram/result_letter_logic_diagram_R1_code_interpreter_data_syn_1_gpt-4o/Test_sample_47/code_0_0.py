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

    # Possible letters
    letters = set("abcdefg")

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    for letter in letters:
        if all(grid[i][6-i] == "" or grid[i][6-i] == letter for i in range(7)):
            # Fill the minor diagonal with this letter
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == "":
                # Find the missing letter for this row
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()