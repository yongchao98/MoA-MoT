# Initial grid setup
grid = [
    ["", "", "g", "d", "", "", ""],
    ["a", "g", "", "c", "", "e", ""],
    ["g", "", "c", "f", "", "", ""],
    ["d", "c", "f", "", "", "a", ""],
    ["", "", "", "", "a", "g", "d"],
    ["f", "e", "", "", "g", "d", ""],
    ["", "b", "a", "g", "", "c", ""]
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    # Check minor diagonal
    if row + col == 6:
        for r in range(7):
            if grid[r][6-r] != "" and grid[r][6-r] != letter:
                return False
    return True

# Backtracking function to fill the grid
def solve(grid):
    for r in range(7):
        for c in range(7):
            if grid[r][c] == "":
                for letter in "abcdefg":
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        if solve(grid):
                            return True
                        grid[r][c] = ""
                return False
    return True

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))