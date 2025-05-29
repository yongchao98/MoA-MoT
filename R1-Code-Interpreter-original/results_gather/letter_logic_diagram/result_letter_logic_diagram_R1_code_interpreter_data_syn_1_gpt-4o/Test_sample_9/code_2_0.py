def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
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

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ""

    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def main():
    grid = [
        ["", "c", "", "", "b", "", ""],
        ["c", "f", "", "", "a", "", ""],
        ["", "", "b", "a", "e", "", "c"],
        ["", "b", "", "", "", "", "f"],
        ["", "", "e", "g", "c", "f", ""],
        ["", "", "g", "", "", "d", ""],
        ["e", "", "", "", "", "b", ""]
    ]

    for letter in 'abcdefg':
        # Try each letter for the minor diagonal
        fill_minor_diagonal(grid, letter)
        if solve(grid):
            for row in grid:
                print(",".join(row))
            return

    print("No solution found")

main()