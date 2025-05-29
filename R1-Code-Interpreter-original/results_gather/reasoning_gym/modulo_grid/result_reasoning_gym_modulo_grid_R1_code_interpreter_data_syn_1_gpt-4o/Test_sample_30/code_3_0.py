# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the pattern
for start_col in range(20):  # Start from each column
    for i in range(7):  # 7 because the diagonal length is 7
        row = 6 - i
        col = start_col + i
        if col < 20 and row >= 0:
            grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))