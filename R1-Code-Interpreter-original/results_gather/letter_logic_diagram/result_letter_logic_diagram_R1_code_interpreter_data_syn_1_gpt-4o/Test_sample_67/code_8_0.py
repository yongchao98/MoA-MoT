def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != "":
        return solve(grid, row, col + 1)
    
    for letter in "abcdefg":
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ""
    return False

def fill_grid_with_diagonal(diagonal_letter):
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Use backtracking to fill the rest of the grid
    if solve(grid):
        return True
    return False

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

# Try each letter for the minor diagonal
for letter in "abcdefg":
    # Make a copy of the grid to try with the current diagonal letter
    grid_copy = [row[:] for row in grid]
    if fill_grid_with_diagonal(letter):
        break

# Print the filled grid
for row in grid:
    print(",".join(row))