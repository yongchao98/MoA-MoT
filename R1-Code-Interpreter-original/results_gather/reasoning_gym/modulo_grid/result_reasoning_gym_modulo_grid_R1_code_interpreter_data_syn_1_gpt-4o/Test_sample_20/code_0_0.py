# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    for j in range(20):
        if i % 6 == 2 and j % 10 == 2:
            grid[i][j] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))