# Define the size of the grid
n = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(n)] for _ in range(n)]

# Fill the grid based on the identified pattern
for i in range(n):
    for j in range(n):
        if (j - i) % 11 == 0:
            grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))