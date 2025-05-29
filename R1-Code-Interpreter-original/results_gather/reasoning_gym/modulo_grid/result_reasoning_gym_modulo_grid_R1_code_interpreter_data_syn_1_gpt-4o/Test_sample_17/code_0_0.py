# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for '✅'
pattern = [2, 3, 0, 1]

# Fill the grid according to the pattern
for i in range(20):
    grid[i][pattern[i % 4]] = '✅'

# Convert the grid to a string format
completed_grid = '\n'.join([''.join(row) for row in grid])

# Print the completed grid
print(completed_grid)