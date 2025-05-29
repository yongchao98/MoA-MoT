# Initialize the grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# List of given positions of ✅
positions = [
    (1, 18), (2, 9), (2, 19), (3, 6), (6, 3), (6, 12), (7, 14),
    (9, 2), (11, 18), (13, 6), (14, 7), (14, 17), (17, 14), (18, 1),
    (18, 10), (19, 2)
]

# Fill the given positions with ✅
for r, c in positions:
    grid[r][c] = '✅'

# Define the pattern of differences
differences = [17, 7, 3, -3, -7, -17]

# Fill the grid based on the pattern
for r in range(20):
    for diff in differences:
        c = r + diff
        if 0 <= c < 20:
            grid[r][c] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))