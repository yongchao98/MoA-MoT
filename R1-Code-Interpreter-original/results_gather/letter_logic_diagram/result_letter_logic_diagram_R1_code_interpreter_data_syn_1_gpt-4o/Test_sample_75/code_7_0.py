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
    for letter in "abcdefg":
        # Check if the letter can be placed on the minor diagonal
        if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            
            # Fill the rest of the grid
            for r in range(7):
                for c in range(7):
                    if grid[r][c] == "":
                        for l in "abcdefg":
                            if is_valid(grid, r, c, l):
                                grid[r][c] = l
                                break
            break

# Initial grid setup
grid = [
    ["", "", "g", "d", "", "", ""],
    ["a", "g", "", "c", "", "e", ""],
    ["g", "", "c", "f", "", "", ""],
    ["d", "c", "f", "", "", "a", ""],
    ["", "", "", "", "a", "g", "d"],
    ["f", "e", "", "", "g", "d", ""],
    ["", "b", "a", "g", "", "c", ""]
]

solve(grid)

# Print the filled grid
for row in grid:
    print(",".join(row))