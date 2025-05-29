from collections import Counter

# Initial grid setup
grid = [
    ["", "d", "g", "f", "a", "b", ""],
    ["d", "", "f", "a", "", "c", "e"],
    ["g", "f", "", "b", "c", "e", ""],
    ["f", "a", "", "c", "e", "d", ""],
    ["a", "b", "", "", "d", "g", "f"],
    ["", "c", "", "d", "", "", "a"],
    ["", "e", "d", "g", "f", "a", "b"]
]

# Function to find missing letters in a list
def find_missing_letters(lst):
    all_letters = set("abcdefg")
    return list(all_letters - set(lst))

# Fill the grid
for i in range(7):
    # Fill rows
    missing_in_row = find_missing_letters(grid[i])
    for j in range(7):
        if grid[i][j] == "":
            grid[i][j] = missing_in_row.pop()

# Check columns and adjust if necessary
for j in range(7):
    column = [grid[i][j] for i in range(7)]
    missing_in_col = find_missing_letters(column)
    for i in range(7):
        if grid[i][j] == "":
            grid[i][j] = missing_in_col.pop()

# Determine the letter for the minor diagonal
minor_diagonal = [grid[i][6-i] for i in range(7)]
counter = Counter(minor_diagonal)
most_common_letter, _ = counter.most_common(1)[0]

# Fill the minor diagonal with the most common letter
for i in range(7):
    grid[i][6-i] = most_common_letter

# Print the completed grid
for row in grid:
    print(",".join(row))