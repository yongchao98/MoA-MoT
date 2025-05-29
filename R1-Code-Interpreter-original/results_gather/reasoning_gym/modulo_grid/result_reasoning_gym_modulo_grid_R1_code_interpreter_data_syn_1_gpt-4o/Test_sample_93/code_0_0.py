# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# List of positions for '✅'
positions = [
    (1, 6), (2, 3), (3, 2), (4, 11), (5, 5), (6, 1), (7, 9), (8, 15),
    (9, 7), (10, 13), (11, 12), (12, 4), (13, 10), (14, 18), (15, 14),
    (16, 8), (17, 17), (18, 16), (19, 19), (5, 20)
]

# Place '✅' in the specified positions
for row, col in positions:
    grid[row-1][col-1] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))