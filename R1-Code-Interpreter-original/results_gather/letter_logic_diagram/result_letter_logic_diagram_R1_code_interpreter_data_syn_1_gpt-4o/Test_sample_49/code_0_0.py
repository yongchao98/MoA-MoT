from itertools import permutations

# Initial grid setup
grid = [
    ["", "", "", "", "d", "", "e"],
    ["c", "f", "", "", "", "g", ""],
    ["f", "b", "", "", "", "a", "c"],
    ["", "d", "", "e", "", "c", ""],
    ["", "g", "", "a", "", "", "b"],
    ["", "", "a", "c", "", "", ""],
    ["", "", "", "", "b", "d", ""]
]

# Function to check if a grid is valid
def is_valid(grid):
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
            return False
    # Check minor diagonal
    minor_diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(minor_diagonal)) != 1:
        return False
    return True

# Try to fill the grid
def fill_grid(grid):
    letters = set("abcdefg")
    # Try all permutations for the minor diagonal
    for letter in letters:
        # Fill the minor diagonal
        for i in range(7):
            grid[i][6-i] = letter
        # Try to fill the rest of the grid
        for perm in permutations(letters):
            for i in range(7):
                if grid[i][i] == "":
                    grid[i][i] = perm[i]
            if is_valid(grid):
                return grid
    return None

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(",".join(row))
else:
    print("No solution found")