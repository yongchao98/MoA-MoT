def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid):
    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, diag_letter) for r in range(7)):
            # Create a copy of the grid to try this configuration
            new_grid = [row[:] for row in grid]
            for r in range(7):
                new_grid[r][6-r] = diag_letter
            if fill_grid(new_grid, 0, 0):
                return new_grid
    return None

def fill_grid(grid, row, col):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != "":
        return fill_grid(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if can_place(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ""
    return False

# Initial grid setup
grid = [
    ["", "", "", "b", "", "e", ""],
    ["f", "", "", "d", "", "c", "a"],
    ["", "", "", "e", "", "a", "f"],
    ["b", "d", "", "", "", "f", ""],
    ["d", "e", "", "", "", "", ""],
    ["", "", "", "f", "g", "", ""],
    ["c", "", "f", "", "b", "", ""]
]

# Solve the puzzle
solution = solve(grid)

# Print the solution
if solution:
    for row in solution:
        print(",".join(row))
else:
    print("No solution found")