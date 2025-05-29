# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the pattern with ✅
for row in range(3, 20, 2):  # Start from the 4th row (index 3) and step by 2
    grid[row][3] = '✅'  # 4th column (index 3)
    grid[row][15] = '✅'  # 16th column (index 15)

# Print the completed grid
for row in grid:
    print(''.join(row))