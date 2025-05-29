def is_valid(grid, row, col, char):
    # Check if char is not in the current row and column
    for i in range(7):
        if grid[row][i] == char or grid[i][col] == char:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for char in "abcdefg":
                    if is_valid(grid, row, col, char):
                        grid[row][col] = char
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = ""
                return False
    return True

# Initial grid setup
grid = [
    ["", "b", "", "c", "d", "f", "a"],
    ["", "", "c", "", "", "", "g"],
    ["", "c", "", "", "", "g", ""],
    ["c", "d", "", "a", "", "", "e"],
    ["d", "", "a", "", "b", "e", ""],
    ["f", "a", "", "", "", "", ""],
    ["a", "", "b", "e", "c", "", "f"]
]

# Choose a letter for the minor diagonal
minor_diagonal_letter = "g"
for i in range(7):
    grid[i][6-i] = minor_diagonal_letter

# Solve the puzzle
solve(grid, minor_diagonal_letter)

# Print the completed grid
for row in grid:
    print(",".join(row))