import numpy as np

# Initial grid setup
grid = [
    ["", "a", "g", "", "", "", ""],
    ["a", "g", "", "f", "b", "", "d"],
    ["", "c", "", "", "", "", ""],
    ["", "f", "b", "", "d", "", ""],
    ["", "", "", "", "", "g", ""],
    ["", "e", "", "", "g", "c", ""],
    ["", "", "a", "g", "c", "", "b"]
]

# Function to check if a letter can be placed in a specific position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Determine the letter for the minor diagonal
# We will try each letter from 'a' to 'g' and see if it can fit
letters = set("abcdefg")
minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

# Find a suitable letter for the minor diagonal
minor_diagonal_letter = None
for letter in letters:
    if all(can_place(grid, r, c, letter) for r, c in minor_diagonal_positions):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal with the chosen letter
for r, c in minor_diagonal_positions:
    grid[r][c] = minor_diagonal_letter

# Function to fill the grid
def fill_grid(grid):
    for r in range(7):
        for c in range(7):
            if grid[r][c] == "":
                for letter in letters:
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        break

# Fill the grid
fill_grid(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))