# Fill the entire grid with checkmarks
grid_size = 20
grid = [['✅' for _ in range(grid_size)] for _ in range(grid_size)]

# Print the completed grid
for row in grid:
    print(''.join(row))