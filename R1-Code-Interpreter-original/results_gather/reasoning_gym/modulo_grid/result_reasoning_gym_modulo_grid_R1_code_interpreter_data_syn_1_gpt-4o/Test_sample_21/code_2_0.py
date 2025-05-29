# Initialize the grid
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern
for i in range(1, 20, 3):
    if (i // 3) % 2 == 0:
        grid[i][7] = '✅'  # Column 8 (index 7)
    else:
        grid[i][3] = '✅'  # Column 4 (index 3)

# Fill in the missing mark
grid[16][15] = '✅'  # Column 16 (index 15)

# Print the grid
for row in grid:
    print(''.join(row))