# Initialize the grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Define the pattern for "✅" marks
pattern = [(3, 13), (9, 19)]

# Fill in the grid based on the pattern
for i in range(20):
    if i % 4 == 1:  # Rows 2, 6, 10, 14, 18
        for pos in pattern[0]:
            grid[i][pos] = "✅"
    elif i % 4 == 3:  # Rows 4, 8, 12, 16, 20
        for pos in pattern[1]:
            grid[i][pos] = "✅"

# Manually set the first row based on the pattern
for pos in pattern[0]:
    grid[0][pos] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))