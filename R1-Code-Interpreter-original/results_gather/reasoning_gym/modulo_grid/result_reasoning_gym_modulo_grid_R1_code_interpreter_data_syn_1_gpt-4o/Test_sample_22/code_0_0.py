# Initialize a 20x20 grid with "❌"
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(1, 11):
    grid[i-1][20-i] = '✅'  # Top-right to bottom-left diagonal
    grid[20-i][i-1] = '✅'  # Bottom-left to top-right diagonal

# Print the completed grid
for row in grid:
    print(''.join(row))