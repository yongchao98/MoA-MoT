# Define the grid
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place the ✅ according to the identified pattern
for i in range(1, 20, 2):
    grid[i][13] = '✅'

# Replace the question mark with ❌
grid[16][10] = '❌'

# Print the completed grid
for row in grid:
    print(''.join(row))