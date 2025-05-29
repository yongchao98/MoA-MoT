from itertools import permutations

# Initial grid setup
grid = [
    ["", "", "", "", "g", "", ""],
    ["", "", "", "", "", "", ""],
    ["e", "", "g", "c", "", "", "b"],
    ["", "g", "", "d", "f", "", "e"],
    ["g", "c", "", "", "b", "", "a"],
    ["c", "d", "", "", "e", "", "g"],
    ["", "f", "", "e", "", "", ""]
]

# Function to check if a letter can be placed in a specific position
def can_place_letter(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Determine the letter for the minor diagonal
letters = set("abcdefg")
minor_diagonal_letter = None
for letter in letters:
    if all(can_place_letter(grid, i, 6-i, letter) for i in range(7)):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal
for i in range(7):
    grid[i][6-i] = minor_diagonal_letter

# Function to fill the grid using backtracking
def fill_grid(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in letters:
                    if can_place_letter(grid, row, col, letter):
                        grid[row][col] = letter
                        if fill_grid(grid):
                            return True
                        grid[row][col] = ""
                return False
    return True

# Attempt to fill the grid
fill_grid(grid)

# Convert all elements to strings for printing
for row in grid:
    print(",".join(map(str, row)))