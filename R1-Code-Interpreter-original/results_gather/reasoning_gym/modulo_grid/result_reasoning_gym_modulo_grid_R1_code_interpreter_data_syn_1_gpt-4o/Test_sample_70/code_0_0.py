# Define the grid with the known pattern
grid = [
    ['❌'] * 20 for _ in range(20)
]

# Fill in the known '✅' positions
positions = [
    (1, 10), (2, 6), (2, 7), (2, 19),
    (5, 4), (5, 17), (7, 10),
    (10, 2), (10, 11), (10, 15),
    (11, 4), (11, 17), (13, 10),
    (14, 6), (14, 7), (14, 19),
    (17, 4), (17, 17), (19, 10)
]

for r, c in positions:
    grid[r][c] = '✅'

# Fill in the missing '❔' with the correct value
grid[18][18] = '❌'

# Print the completed grid
for row in grid:
    print(''.join(row))