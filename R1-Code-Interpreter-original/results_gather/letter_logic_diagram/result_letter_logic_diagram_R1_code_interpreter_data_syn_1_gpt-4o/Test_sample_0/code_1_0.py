def solve_puzzle():
    # Initial grid setup
    grid = [
        ["", "", "e", "", "", "", ""],
        ["g", "e", "f", "", "", "", "c"],
        ["", "", "", "a", "", "", ""],
        ["", "", "", "", "", "", "e"],
        ["", "a", "d", "c", "", "", ""],
        ["", "d", "c", "g", "", "", "b"],
        ["d", "c", "g", "", "f", "b", ""]
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    possible_letters = set("abcdefg")
    for i, j in zip(range(7), range(6, -1, -1)):
        if grid[i][j] != "":
            possible_letters.intersection_update(grid[i][j])

    # Choose the letter for the minor diagonal
    diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i, j in zip(range(7), range(6, -1, -1)):
        grid[i][j] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set("abcdefg") - set(grid[i])
        for j in range(7):
            if grid[i][j] == "":
                # Find a letter that can be placed in this position
                col_letters = set(grid[k][j] for k in range(7))
                possible_letters = row_letters - col_letters
                if possible_letters:
                    chosen_letter = possible_letters.pop()
                    grid[i][j] = chosen_letter
                    row_letters.remove(chosen_letter)
                else:
                    # If no possible letters, print an error message
                    print(f"Error: No available letters for cell ({i}, {j})")
                    return

    # Print the final grid
    for row in grid:
        print(",".join(row))

solve_puzzle()