from itertools import permutations

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

# Function to solve the grid
def solve(grid):
    # Determine the letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Try to fill the grid with this letter on the diagonal
            if fill_grid(grid, letter):
                return grid
    return None

# Function to fill the grid with backtracking
def fill_grid(grid, diag_letter):
    # Fill the minor diagonal
    for r in range(7):
        grid[r][6-r] = diag_letter
    
    # Try to fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == "":
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        if fill_grid(grid, diag_letter):
                            return True
                        grid[r][c] = ""
                return False
    return True

# Solve the puzzle
solution = solve(grid)

# Print the solution
if solution:
    for row in solution:
        print(",".join(row))
else:
    print("No solution found")