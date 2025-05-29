def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ["", "", "d", "g", "", "b", "f"],
        ["", "d", "", "", "", "f", ""],
        ["d", "g", "", "", "f", "", ""],
        ["g", "", "b", "f", "e", "c", "d"],
        ["a", "b", "f", "", "c", "d", ""],
        ["", "", "", "c", "d", "g", ""],
        ["f", "", "c", "d", "g", "a", "b"]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    possible_letters = set("abcdefg")
    for i in range(7):
        for j in range(7):
            if grid[i][j] != "":
                possible_letters.discard(grid[i][j])

    # Choose a letter for the minor diagonal
    diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == "":
                # Find the missing letter for this row
                row_letters = set(grid[i])
                missing_letter = (set("abcdefg") - row_letters).pop()
                grid[i][j] = missing_letter

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()