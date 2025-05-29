def solve_puzzle():
    # Initial grid setup
    grid = [
        ["", "", "c", "f", "", "", ""],
        ["e", "c", "f", "", "", "", ""],
        ["", "f", "g", "b", "a", "", "e"],
        ["f", "g", "b", "a", "", "", "c"],
        ["g", "b", "a", "d", "e", "", "f"],
        ["", "", "", "", "c", "", ""],
        ["a", "", "e", "", "", "", ""]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is most suitable for the diagonal
    # We will choose 'd' for the diagonal as it is not present in any of the diagonal positions
    diagonal_letter = 'd'
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Fill the rest of the grid
    letters = set("abcdefg")
    for i in range(7):
        for j in range(7):
            if grid[i][j] == "":
                # Determine which letters are already used in the row and column
                used_letters = set(grid[i]) | {grid[k][j] for k in range(7)}
                # Find the missing letter
                missing_letter = (letters - used_letters).pop()
                grid[i][j] = missing_letter

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()