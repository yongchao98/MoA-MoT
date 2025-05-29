# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Define the pattern for "✅"
pattern_1 = [3, 10, 17]  # Columns for rows 1, 8, 14, 19 (0-based index)
pattern_2 = [5, 12, 19]  # Columns for rows 6, 13, 18 (0-based index)

# Fill the grid with "✅" based on the pattern
for i in range(20):
    if i % 5 == 0 or i % 5 == 3:  # Rows 1, 8, 14, 19 (0-based index)
        for col in pattern_1:
            grid[i][col] = "✅"
    elif i % 5 == 1:  # Rows 6, 13, 18 (0-based index)
        for col in pattern_2:
            grid[i][col] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))