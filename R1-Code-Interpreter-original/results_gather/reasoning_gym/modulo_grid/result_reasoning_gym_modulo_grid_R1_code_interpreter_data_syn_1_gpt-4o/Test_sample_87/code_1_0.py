# Define the grid size
rows, cols = 20, 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Define the pattern for '✅'
pattern = [
    [6, 15],  # Row 0 pattern
    [6],      # Row 1 pattern
    [2, 4, 6, 8, 12, 14, 16],  # Row 2 pattern
    [6]       # Row 3 pattern
]

# Fill the grid based on the observed pattern
for r in range(rows):
    for c in pattern[r % 4]:
        grid[r][c] = '✅'

# Convert the grid to a string for display
grid_str = '\n'.join([''.join(row) for row in grid])
print(grid_str)