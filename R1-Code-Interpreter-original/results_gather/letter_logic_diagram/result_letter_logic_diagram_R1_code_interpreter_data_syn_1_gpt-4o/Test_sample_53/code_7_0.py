def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid):
    # Try each letter for the minor diagonal
    for diag_letter in "abcdefg":
        # Check if the letter can be placed on the minor diagonal
        if all(is_valid(grid, i, 6-i, diag_letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = diag_letter
            # Try to fill the rest of the grid
            if fill_grid(grid, 0, 0):
                return grid
            # Reset the minor diagonal if it doesn't work
            for i in range(7):
                grid[i][6-i] = ""
    return None

def fill_grid(grid, row, col):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != "":
        return fill_grid(grid, row, col + 1)
    for letter in "abcdefg":
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ""
    return False

# Initial grid setup
grid = [
    ["", "", "f", "", "", "", ""],
    ["", "", "d", "e", "g", "", ""],
    ["", "", "e", "", "", "", "a"],
    ["", "", "", "c", "b", "a", ""],
    ["e", "", "", "b", "", "", "f"],
    ["", "c", "", "a", "f", "", ""],
    ["", "", "", "f", "d", "e", "g"]
]

# Solve the grid
solution = solve(grid)

# Print the solution
if solution:
    for row in solution:
        print(",".join(row))
else:
    print("No solution found")