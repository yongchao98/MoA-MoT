def fill_grid(grid):
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

    # Fill the minor diagonal with the letter 'e'
    for i in range(7):
        grid[i][6-i] = 'e'

    # Fill the rest of the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in "abcdefg":
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        break

    return grid

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

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(",".join(row))