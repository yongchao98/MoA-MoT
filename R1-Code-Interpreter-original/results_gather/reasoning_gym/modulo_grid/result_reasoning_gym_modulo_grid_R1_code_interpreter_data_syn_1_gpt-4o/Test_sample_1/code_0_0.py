# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    grid[i][i] = "✅"

# Convert the grid to the required format
completed_grid = "\n".join(["".join(row) for row in grid])

print(completed_grid)