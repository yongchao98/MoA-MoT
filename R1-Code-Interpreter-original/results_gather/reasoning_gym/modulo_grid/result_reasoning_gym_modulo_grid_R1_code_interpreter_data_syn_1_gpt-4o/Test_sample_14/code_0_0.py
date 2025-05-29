# Initialize a 20x20 grid with "❌"
grid = [["❌" for _ in range(20)] for _ in range(20)]

# Fill in the "✅" based on the identified pattern
for i in range(20):
    if i < 6:
        grid[i][5 - i] = "✅"
    else:
        grid[i][i - 6] = "✅"

# Print the completed grid
for row in grid:
    print("".join(row))