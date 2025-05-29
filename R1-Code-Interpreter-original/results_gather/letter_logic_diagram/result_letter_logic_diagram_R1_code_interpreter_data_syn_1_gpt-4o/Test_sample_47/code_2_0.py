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

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == "":
                    for letter in letters:
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            if solve(grid):
                                return True
                            grid[i][j] = ""
                    return False
        return True

    solve(grid)

    # Print the completed grid
    for row in grid:
        print(",".join(row))

solve_puzzle()