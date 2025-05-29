# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the checkmarks based on the identified pattern
for row in range(20):
    if row % 3 == 2:  # Pattern observed in Row 3
        for col in range(3, 20, 4):
            grid[row][col] = '✅'

# Fill in the question mark position
grid[17][0] = '✅'  # Row 18, Column 1 (0-indexed)

# Print the completed grid
for row in grid:
    print(''.join(row))