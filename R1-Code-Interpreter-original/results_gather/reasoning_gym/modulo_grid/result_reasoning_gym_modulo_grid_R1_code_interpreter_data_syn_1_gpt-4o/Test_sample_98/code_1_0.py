# Define the positions of the "✅" marks
positions = [
    (1, 15), (2, 8), (3, 17), (4, 13), (5, 12), (6, 14), (7, 9),
    (8, 3), (8, 20), (9, 7), (10, 11), (11, 16), (12, 10), (13, 5),
    (14, 6), (15, 2), (15, 19), (16, 10), (17, 4), (19, 15), (20, 8)
]

# Create a 20x20 grid initialized with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill the grid with "✅" based on the positions
for row, col in positions:
    grid[row - 1][col - 1] = "✅"

# Identify the missing position
# Check each row to find the missing "✅"
for row in range(20):
    if grid[row].count("✅") < 2:
        missing_col = grid[row].index("❌")
        grid[row][missing_col] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))