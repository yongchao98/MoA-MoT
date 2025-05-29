# Initialize the grid
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the first row with ✅
for i in range(20):
    grid[0][i] = '✅'

# Fill the 5th, 9th, 13th, and 17th rows with alternating pattern
for row in [4, 8, 12, 16]:
    for col in range(20):
        if col % 2 == 1:
            grid[row][col] = '✅'

# Fill the 2nd and 18th columns with ✅
for row in range(20):
    grid[row][1] = '✅'
    grid[row][17] = '✅'

# Fill the 9th column of the 2nd, 6th, 10th, 14th, and 18th rows with ✅
for row in [1, 5, 9, 13, 17]:
    grid[row][8] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))