from itertools import permutations

# Initial grid setup
grid = [
    ["", "c", "", "e", "", "", ""],
    ["c", "", "e", "", "", "", "d"],
    ["a", "", "b", "f", "", "", "c"],
    ["", "", "", "", "", "c", "a"],
    ["b", "", "g", "d", "", "", ""],
    ["", "g", "", "", "", "e", "b"],
    ["", "", "c", "", "e", "", ""]
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
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            
            # Try to fill the rest of the grid
            if fill_remaining(grid):
                return grid
            else:
                # Reset the minor diagonal if it doesn't work
                for i in range(7):
                    grid[i][6-i] = ""
    return None

# Function to fill the remaining grid
def fill_remaining(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        if fill_remaining(grid):
                            return True
                        grid[row][col] = ""
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