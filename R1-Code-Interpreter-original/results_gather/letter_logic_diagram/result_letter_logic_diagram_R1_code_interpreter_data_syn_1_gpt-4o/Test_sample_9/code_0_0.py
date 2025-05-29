from itertools import permutations

# Initial grid setup
grid = [
    ["", "c", "", "", "b", "", ""],
    ["c", "f", "", "", "a", "", ""],
    ["", "", "b", "a", "e", "", "c"],
    ["", "b", "", "", "", "", "f"],
    ["", "", "e", "g", "c", "f", ""],
    ["", "", "g", "", "", "d", ""],
    ["e", "", "", "", "", "b", ""]
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
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            
            # Try to fill the rest of the grid
            if fill_remaining(grid):
                return grid
            
            # Reset the minor diagonal if it doesn't work
            for r in range(7):
                grid[r][6-r] = ""
    return None

# Function to fill the remaining grid
def fill_remaining(grid):
    for r in range(7):
        for c in range(7):
            if grid[r][c] == "":
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        if fill_remaining(grid):
                            return True
                        grid[r][c] = ""
                return False
    return True

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(",".join(row))
else:
    print("No solution found")