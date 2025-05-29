# Initialize the grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with ✅ based on the identified pattern
for row in range(20):
    if row < 10:
        col = 10 - row
    else:
        col = row - 10
    if 0 <= col < 20:
        grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))