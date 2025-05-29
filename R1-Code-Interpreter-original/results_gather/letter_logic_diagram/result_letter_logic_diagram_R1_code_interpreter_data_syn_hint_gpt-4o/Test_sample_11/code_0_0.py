from itertools import permutations

# Initial grid setup
grid = [
    ["", "", "d", "g", "", "b", "f"],
    ["", "d", "", "", "", "f", ""],
    ["d", "g", "", "", "f", "", ""],
    ["g", "", "b", "f", "e", "c", "d"],
    ["a", "b", "f", "", "c", "d", ""],
    ["", "", "", "c", "d", "g", ""],
    ["f", "", "c", "d", "g", "a", "b"]
]

# Function to check if a letter can be placed at grid[row][col]
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to solve the grid
def solve(grid):
    # Try to fill the minor diagonal with the same letter
    for letter in "abcdefg":
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Fill the rest of the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in "abcdefg":
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid):
                            return True
                        grid[row][col] = ""
                return False
    return True

# Solve the grid
solve(grid)

# Print the solved grid
for row in grid:
    print(",".join(row))