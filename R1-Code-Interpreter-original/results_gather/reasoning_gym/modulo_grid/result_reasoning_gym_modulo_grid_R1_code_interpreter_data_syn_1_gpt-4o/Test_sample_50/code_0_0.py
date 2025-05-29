# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' based on the identified pattern
for row in range(1, 20, 2):  # Start from row 2 (index 1) and step by 2
    for col in [2, 8, 14]:  # Columns 3, 9, 15 (0-indexed as 2, 8, 14)
        grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))