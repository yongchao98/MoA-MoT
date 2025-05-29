from itertools import permutations

# Initial grid setup
grid = [
    ["", "", "d", "c", "a", "", ""],
    ["", "d", "c", "a", "", "f", ""],
    ["d", "", "", "", "", "", "g"],
    ["", "", "", "", "", "", "d"],
    ["", "", "", "", "", "", ""],
    ["b", "f", "", "g", "d", "", "a"],
    ["f", "", "g", "", "", "a", "b"]
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
    return True

# Backtracking function to fill the grid
def fill_grid(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != "":
        return fill_grid(grid, row, col + 1)

    # Try each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if can_place(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ""

    return False

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all(can_place(grid, r, 6-r, letter) for r in range(7)):
        # Place the letter on the minor diagonal
        for r in range(7):
            grid[r][6-r] = letter
        break

# Fill the grid
fill_grid(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))