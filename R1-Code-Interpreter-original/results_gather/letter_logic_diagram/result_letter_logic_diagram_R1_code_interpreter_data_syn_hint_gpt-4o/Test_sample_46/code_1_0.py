def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Try to fill the grid using backtracking
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

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

# Initial grid setup
grid = [
    ["", "g", "e", "", "b", "", ""],
    ["g", "e", "f", "b", "", "", ""],
    ["", "", "", "", "", "", "g"],
    ["f", "", "d", "", "", "g", "e"],
    ["", "d", "", "", "", "", ""],
    ["d", "c", "", "", "", "", "b"],
    ["c", "", "g", "", "", "b", ""]
]

# Try each letter for the minor diagonal
for letter in "abcdefg":
    # Make a copy of the grid to try this letter
    test_grid = [row[:] for row in grid]
    fill_minor_diagonal(test_grid, letter)
    if solve(test_grid, letter):
        for row in test_grid:
            print(",".join(row))
        break