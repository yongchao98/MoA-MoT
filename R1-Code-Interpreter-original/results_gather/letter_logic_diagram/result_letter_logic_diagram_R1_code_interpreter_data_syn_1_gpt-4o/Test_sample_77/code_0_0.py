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

# Function to fill the grid
def fill_grid(grid):
    # Determine the letter for the minor diagonal
    # Try each letter from 'a' to 'g' and see if it fits
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            break

    # Fill the rest of the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        break

    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(",".join(row))