def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col] without conflicts
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in "abcdefg":
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = ""
                return False
    return True

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
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set("abcdefg")
    minor_diagonal_letters = set(grid[i][j] for i, j in minor_diagonal_indices if grid[i][j])
    if len(minor_diagonal_letters) == 1:
        minor_diagonal_letter = minor_diagonal_letters.pop()
    else:
        minor_diagonal_letter = list(letters - minor_diagonal_letters)[0]

    # Fill the minor diagonal with the determined letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter

    # Solve the puzzle
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(",".join(row))
    else:
        print("No solution found")

solve_puzzle()