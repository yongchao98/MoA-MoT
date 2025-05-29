# Initialize the grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Define the pattern for "✅"
pattern = [
    (1, 3), (1, 13),
    (3, 9), (3, 19),
    (5, 3), (5, 13),
    (7, 9), (7, 19),
    (9, 3), (9, 13),
    (11, 9), (11, 19),
    (13, 3), (13, 13),
    (15, 9), (15, 19),
    (17, 3), (17, 13),
    (19, 9), (19, 19)
]

# Fill the grid based on the pattern
for row, col in pattern:
    grid[row][col] = "✅"

# Fill the row with the question mark
grid[0][9] = "✅"
grid[0][19] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))