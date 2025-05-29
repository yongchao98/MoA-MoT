from itertools import permutations

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

# Function to check if a letter can be placed on the minor diagonal
def can_place_on_diagonal(letter):
    for i in range(7):
        if grid[i][6-i] not in ("", letter):
            return False
    return True

# Function to fill the grid
def fill_grid(diagonal_letter):
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        col_letters = set(grid[j][i] for j in range(7))
        missing_letters = set("abcdefg") - row_letters
        for j in range(7):
            if grid[i][j] == "":
                for letter in missing_letters:
                    if letter not in col_letters:
                        grid[i][j] = letter
                        col_letters.add(letter)
                        missing_letters.remove(letter)
                        break

# Try each letter for the minor diagonal
for letter in "abcdefg":
    if can_place_on_diagonal(letter):
        fill_grid(letter)
        break

# Print the filled grid
for row in grid:
    print(",".join(row))