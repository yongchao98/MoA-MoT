# Initial grid setup
grid = [
    ["", "b", "", "c", "d", "f", "a"],
    ["", "", "c", "", "", "", "g"],
    ["", "c", "", "", "", "g", ""],
    ["c", "d", "", "a", "", "", "e"],
    ["d", "", "a", "", "b", "e", ""],
    ["f", "a", "", "", "", "", ""],
    ["a", "", "b", "e", "c", "", "f"]
]

# Function to find the missing letter in a list
def find_missing_letters(lst):
    all_letters = set("abcdefg")
    return list(all_letters - set(lst))

# Fill the minor diagonal with the same letter
minor_diagonal_letter = "g"  # Choosing 'g' as it appears frequently
for i in range(7):
    grid[i][6-i] = minor_diagonal_letter

# Fill in the rest of the grid
for i in range(7):
    # Fill rows
    missing_letters = find_missing_letters(grid[i])
    for j in range(7):
        if grid[i][j] == "":
            grid[i][j] = missing_letters.pop()

    # Fill columns
    for j in range(7):
        column = [grid[x][j] for x in range(7)]
        missing_letters = find_missing_letters(column)
        for x in range(7):
            if grid[x][j] == "":
                grid[x][j] = missing_letters.pop()

# Print the completed grid
for row in grid:
    print(",".join(row))