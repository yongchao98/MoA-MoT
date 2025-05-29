# Initialize the grid with empty values
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for each set of 6 rows
patterns = [
    [3],  # Row 1 pattern
    [6],  # Row 2 pattern
    [12], # Row 3 pattern
    [1, 4, 7, 10, 13, 16, 19], # Row 4 pattern
    [3],  # Row 5 pattern
    [6]   # Row 6 pattern
]

# Fill the grid based on the pattern
for i in range(20):
    pattern = patterns[i % 6]
    for col in pattern:
        grid[i][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))