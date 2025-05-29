from itertools import permutations

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

# Function to fill the grid
def fill_grid(grid):
    # Try each letter for the minor diagonal
    for letter in "abcdefg":
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            
            # Fill the rest of the grid
            for r in range(7):
                for c in range(7):
                    if grid[r][c] == "":
                        for l in "abcdefg":
                            if can_place(grid, r, c, l):
                                grid[r][c] = l
                                break
            break

fill_grid(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))