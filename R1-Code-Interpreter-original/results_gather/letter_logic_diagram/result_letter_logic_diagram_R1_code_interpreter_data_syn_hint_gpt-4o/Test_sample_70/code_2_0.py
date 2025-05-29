def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, diagonal_positions, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, diagonal_letter, diagonal_positions, row + 1, 0)
    if (row, col) in diagonal_positions:
        if grid[row][col] == diagonal_letter or grid[row][col] == "":
            grid[row][col] = diagonal_letter
            return solve(grid, diagonal_letter, diagonal_positions, row, col + 1)
        else:
            return False
    if grid[row][col] != "":
        return solve(grid, diagonal_letter, diagonal_positions, row, col + 1)
    for letter in "abcdefg":
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter, diagonal_positions, row, col + 1):
                return True
            grid[row][col] = ""
    return False

def solve_puzzle():
    grid = [
        ["", "", "c", "f", "", "", ""],
        ["e", "c", "f", "", "", "", ""],
        ["", "f", "g", "b", "a", "", "e"],
        ["f", "g", "b", "a", "", "", "c"],
        ["g", "b", "a", "d", "e", "", "f"],
        ["", "", "", "", "c", "", ""],
        ["a", "", "e", "", "", "", ""]
    ]

    # Try different letters for the diagonal
    for diagonal_letter in "abcdefg":
        diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
        if solve(grid, diagonal_letter, diagonal_positions):
            for row in grid:
                print(",".join(row))
            return

    print("No solution found")

solve_puzzle()