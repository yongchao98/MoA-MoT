def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != "":
        return solve(grid, row, col + 1)

    for letter in "abcdefg":
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ""

    return False

# Initial grid setup
grid = [
    ["", "g", "b", "", "", "", ""],
    ["", "", "", "", "c", "", ""],
    ["b", "", "", "", "", "", "a"],
    ["f", "d", "c", "", "", "", ""],
    ["", "c", "e", "a", "", "b", ""],
    ["", "", "", "g", "", "f", "d"],
    ["e", "a", "g", "", "f", "", ""]
]

# Determine the letter for the minor diagonal
minor_diagonal_letter = None
for letter in "abcdefg":
    if all(grid[r][c] in ("", letter) for r, c in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal with the chosen letter
for r, c in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
    grid[r][c] = minor_diagonal_letter

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))