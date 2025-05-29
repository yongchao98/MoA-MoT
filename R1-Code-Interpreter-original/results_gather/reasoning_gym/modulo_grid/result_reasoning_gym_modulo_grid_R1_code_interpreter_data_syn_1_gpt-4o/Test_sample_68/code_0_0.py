# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill the grid with "✅" where the sum of indices equals 19
for i in range(20):
    for j in range(20):
        if i + j == 19:
            grid[i][j] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))