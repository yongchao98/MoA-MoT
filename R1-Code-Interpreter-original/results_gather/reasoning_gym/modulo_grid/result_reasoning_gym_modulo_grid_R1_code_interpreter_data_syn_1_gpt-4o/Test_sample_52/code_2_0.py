# Initialize the grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Define the pattern for "✅" marks
pattern = [(3, 13), (9, 19)]

# Fill in the grid based on the pattern
for i in range(20):
    if i % 4 == 0:  # Rows 1, 5, 9, 13, 17
        for pos in pattern[0]:
            grid[i][pos] = "✅"
    elif i % 4 == 2:  # Rows 3, 7, 11, 15, 19
        for pos in pattern[1]:
            grid[i][pos] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))