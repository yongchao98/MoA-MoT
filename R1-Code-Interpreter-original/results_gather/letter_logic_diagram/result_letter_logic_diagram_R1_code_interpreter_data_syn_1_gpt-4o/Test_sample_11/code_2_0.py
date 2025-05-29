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
    diagonal_letter = None
    for letter in "abcdefg":
        can_place = True
        for i in range(7):
            if letter in grid[i] or letter in [grid[j][6-i] for j in range(7)]:
                can_place = False
                break
        if can_place:
            diagonal_letter = letter
            break

    if diagonal_letter is None:
        raise ValueError("No suitable letter found for the minor diagonal")

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == "":
                # Find the missing letter for this row
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (set("abcdefg") - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()