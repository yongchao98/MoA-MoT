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

# Identify the pattern and fill in the '❔' mark
# The pattern seems to be a diagonal or some arithmetic sequence
# Let's check the differences between row and column indices
differences = [(row - col) for row, col in check_marks]

# Print the differences to identify the pattern
print(differences)

# Based on the pattern, fill in the '❔' mark
# Assuming the pattern is consistent, find the missing '✅' position
for i in range(grid_size):
    for j in range(grid_size):
        if grid[i][j] == '❌':
            if (i - j) in differences:
                grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))