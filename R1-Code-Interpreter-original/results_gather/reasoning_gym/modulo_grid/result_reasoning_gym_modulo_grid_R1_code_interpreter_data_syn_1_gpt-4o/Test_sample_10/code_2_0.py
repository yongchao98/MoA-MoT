# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill the grid with "✅" based on the identified pattern
for row in range(1, 20, 4):
    for col in range(3, 20, 4):
        grid[row][col] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))