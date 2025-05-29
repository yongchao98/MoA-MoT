# Define the size of the grid
n = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(n)] for _ in range(n)]

# Fill the grid based on the identified pattern
for i in range(n):
    grid[i][i % 11] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))