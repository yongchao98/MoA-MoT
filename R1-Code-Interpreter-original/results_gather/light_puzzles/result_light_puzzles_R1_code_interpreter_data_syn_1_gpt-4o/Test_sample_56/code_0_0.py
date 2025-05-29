# Initial grid setup
grid = [
    [0, 0],
    [1, 1]
]

# Button A controls these positions
button_a_positions = [(1, 0), (0, 0)]  # (row, column) format

# Toggle function
def toggle_light(grid, row, col):
    grid[row][col] = 1 - grid[row][col]

# Apply Button A
for pos in button_a_positions:
    toggle_light(grid, pos[0], pos[1])

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)