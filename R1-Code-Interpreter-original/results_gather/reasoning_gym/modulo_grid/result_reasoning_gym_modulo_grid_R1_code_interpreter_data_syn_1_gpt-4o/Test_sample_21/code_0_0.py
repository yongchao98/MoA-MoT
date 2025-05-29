# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for ✅ marks
pattern_positions = [(1, 8), (4, 4), (6, 8), (8, 4), (11, 8), (13, 4), (16, 8), (18, 4)]

# Place the ✅ marks according to the pattern
for row, col in pattern_positions:
    grid[row][col] = '✅'

# Place the question mark
grid[16][15] = '❔'

# Convert the grid to a string representation
completed_grid = '\n'.join([''.join(row) for row in grid])

print(completed_grid)