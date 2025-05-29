# Initial grid setup
grid = [
    ["", "c", "b", "", "", "d", ""],
    ["c", "", "g", "e", "d", "f", "a"],
    ["b", "", "", "", "", "", "c"],
    ["", "", "", "", "a", "", ""],
    ["", "", "f", "", "", "", ""],
    ["", "f", "", "c", "", "g", ""],
    ["f", "", "", "b", "g", "", ""]
]

# Function to check if a letter can be placed in a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to fill the minor diagonal with a valid letter
def fill_minor_diagonal(grid):
    # Try each letter for the minor diagonal
    for letter in "abcdefg":
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            for r in range(7):
                grid[r][6-r] = letter
            return True
    return False

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

# Fill the minor diagonal
fill_minor_diagonal(grid)

# Solve the grid
solve(grid)

# Print the final grid
for row in grid:
    print(",".join(row))