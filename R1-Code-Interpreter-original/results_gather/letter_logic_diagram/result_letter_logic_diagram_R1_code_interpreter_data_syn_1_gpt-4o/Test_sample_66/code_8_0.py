def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != "":
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ""

        return False

    if backtrack(0, 0):
        return grid
    else:
        return None

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

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    # Make a deep copy of the grid to try this letter
    grid_copy = [row[:] for row in grid]
    result = solve(grid_copy, letter)
    if result:
        for row in result:
            print(",".join(row))
        break
else:
    print("No solution found")