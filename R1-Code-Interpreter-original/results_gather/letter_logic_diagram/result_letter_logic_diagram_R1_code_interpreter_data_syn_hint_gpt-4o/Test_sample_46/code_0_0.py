from itertools import permutations

# Initial grid setup
grid = [
    ["", "g", "e", "", "b", "", ""],
    ["g", "e", "f", "b", "", "", ""],
    ["", "", "", "", "", "", "g"],
    ["f", "", "d", "", "", "g", "e"],
    ["", "d", "", "", "", "", ""],
    ["d", "c", "", "", "", "", "b"],
    ["c", "", "g", "", "", "b", ""]
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

# Try to fill the minor diagonal with each letter and check for a valid solution
for letter in "abcdefg":
    # Place the letter on the minor diagonal
    for i in range(7):
        grid[i][6-i] = letter
    
    # Try to fill the rest of the grid
    valid = True
    for row in range(7):
        for col in range(7):
            if grid[row][col] == "":
                # Find a valid letter for this position
                for l in "abcdefg":
                    if can_place(grid, row, col, l):
                        grid[row][col] = l
                        break
                else:
                    valid = False
                    break
        if not valid:
            break
    
    if valid:
        break

# Print the completed grid
for row in grid:
    print(",".join(row))