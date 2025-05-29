# Define the grid size
grid_size = 20

# Initialize the grid with crosses
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Define the pattern based on the observed positions
pattern = [
    [8],
    [5, 11, 19],
    [12],
    [3, 10, 17],
    [10],
    [7, 14],
    [],
    [2, 9, 16],
    [5, 19],
    [6, 13, 20],
    [3, 17],
    [4, 11, 18],
    [7],
    [],
    [],
    [8],
    [5, 11, 19],
    [12],
    [3, 10, 17],
    []
]

# Fill the grid with the pattern
for row in range(grid_size):
    for col in pattern[row]:
        grid[row][col - 1] = '✅'  # Adjust for 0-based index

# Print the completed grid
for row in grid:
    print(''.join(row))