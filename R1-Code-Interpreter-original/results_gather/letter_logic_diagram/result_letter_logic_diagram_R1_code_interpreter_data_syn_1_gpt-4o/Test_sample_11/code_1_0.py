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

    # Try each letter for the minor diagonal
    for letter in "abcdefg":
        # Check if the letter can be placed on the diagonal
        can_place = True
        for i in range(7):
            if letter in grid[i] or letter in [grid[j][6-i] for j in range(7)]:
                can_place = False
                break
        if not can_place:
            continue

        # Fill the minor diagonal with the chosen letter
        for i in range(7):
            grid[i][6-i] = letter

        # Fill the rest of the grid
        for i in range(7):
            for j in range(7):
                if grid[i][j] == "":
                    # Find the missing letter for this row and column
                    row_letters = set(grid[i])
                    col_letters = set(grid[k][j] for k in range(7))
                    missing_letter = (set("abcdefg") - row_letters - col_letters).pop()
                    grid[i][j] = missing_letter

        # Check if the grid is valid
        valid = True
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(grid[j][i] for j in range(7))) != 7:
                valid = False
                break

        if valid:
            # Print the completed grid
            for row in grid:
                print(",".join(row))
            return

solve_puzzle()