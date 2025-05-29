# Define the grid size
grid_size = 20

# Define the positions of the '✅' marks
check_marks = [
    (2, 8), (3, 11), (4, 12), (5, 5), (6, 4),
    (7, 13), (8, 2), (9, 10), (10, 9), (11, 3),
    (12, 4), (13, 7), (14, 19)
]

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Place the '✅' marks in the grid
for row, col in check_marks:
    grid[row - 1][col - 1] = '✅'

# Manually inspect the grid to identify the pattern
# Based on the differences, fill in the '❔' mark
# Assuming the pattern is consistent, find the missing '✅' position

# Print the completed grid
for row in grid:
    print(''.join(row))