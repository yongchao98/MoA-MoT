from collections import Counter

# Initial grid setup
grid = [
    ["", "g", "b", "", "", "", ""],
    ["", "", "", "", "c", "", ""],
    ["b", "", "", "", "", "", "a"],
    ["f", "d", "c", "", "", "", ""],
    ["", "c", "e", "a", "", "b", ""],
    ["", "", "", "g", "", "f", "d"],
    ["e", "a", "g", "", "f", "", ""]
]

# Function to check if a letter can be placed on the minor diagonal
def can_place_on_minor_diagonal(letter):
    positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for r, c in positions:
        if grid[r][c] != "" and grid[r][c] != letter:
            return False
    return True

# Find a suitable letter for the minor diagonal
letters = "abcdefg"
minor_diagonal_letter = None
for letter in letters:
    if can_place_on_minor_diagonal(letter):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal with the chosen letter
positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
for r, c in positions:
    grid[r][c] = minor_diagonal_letter

# Function to fill the grid ensuring each row and column has all letters
def fill_grid():
    for r in range(7):
        row_letters = set(grid[r])
        missing_letters = set(letters) - row_letters
        for c in range(7):
            if grid[r][c] == "":
                for letter in missing_letters:
                    # Check if the letter can be placed in the column
                    if all(grid[i][c] != letter for i in range(7)):
                        grid[r][c] = letter
                        missing_letters.remove(letter)
                        break

fill_grid()

# Print the filled grid
for row in grid:
    print(",".join(row))