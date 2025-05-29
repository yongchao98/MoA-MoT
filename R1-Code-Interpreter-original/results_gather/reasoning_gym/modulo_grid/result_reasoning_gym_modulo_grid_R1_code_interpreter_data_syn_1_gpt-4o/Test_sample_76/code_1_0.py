# Initialize the grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the column sets for the pattern
columns_set_1 = [3, 8, 13, 18]  # 0-based index for columns 4, 9, 14, 19
columns_set_2 = [2, 7, 12, 17]  # 0-based index for columns 3, 8, 13, 18

# Fill the grid based on the identified pattern
for row in range(20):
    if row % 4 == 1:
        for col in columns_set_1:
            grid[row][col] = '✅'
    elif row % 4 == 3:
        for col in columns_set_2:
            grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))