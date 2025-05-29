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

def fill_grid():
    grid = [
        ["", "b", "", "", "", "", ""],
        ["b", "a", "g", "d", "", "", "f"],
        ["", "", "d", "", "", "", "b"],
        ["g", "", "", "e", "", "", ""],
        ["", "", "", "", "b", "", ""],
        ["c", "", "f", "", "a", "", "d"],
        ["e", "f", "", "a", "g", "", ""]
    ]

    # Determine the letter for the minor diagonal
    for letter in "abcdefg":
        if all(is_valid(grid, i, 6-i, letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Solve the grid
    if solve(grid):
        for row in grid:
            print(",".join(row))
    else:
        print("No solution found")

fill_grid()