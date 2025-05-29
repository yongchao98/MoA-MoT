def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        # If we've filled all rows, the grid is complete
        return True
    if col == 7:
        # Move to the next row
        return solve(grid, row + 1, 0)
    if grid[row][col] != "":
        # Skip pre-filled cells
        return solve(grid, row, col + 1)

    # Try placing each letter in the current cell
    for letter in "abcdefg":
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ""  # Backtrack

    return False

def solve_puzzle():
    # Initial grid setup
    grid = [
        ["", "", "e", "", "", "", ""],
        ["g", "e", "f", "", "", "", "c"],
        ["", "", "", "a", "", "", ""],
        ["", "", "", "", "", "", "e"],
        ["", "a", "d", "c", "", "", ""],
        ["", "d", "c", "g", "", "", "b"],
        ["d", "c", "g", "", "f", "b", ""]
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set("abcdefg")
    for i, j in zip(range(7), range(6, -1, -1)):
        if grid[i][j] != "":
            possible_letters.intersection_update(grid[i][j])

    diagonal_letter = possible_letters.pop()

    for i, j in zip(range(7), range(6, -1, -1)):
        grid[i][j] = diagonal_letter

    if solve(grid):
        for row in grid:
            print(",".join(row))
    else:
        print("No solution found")

solve_puzzle()