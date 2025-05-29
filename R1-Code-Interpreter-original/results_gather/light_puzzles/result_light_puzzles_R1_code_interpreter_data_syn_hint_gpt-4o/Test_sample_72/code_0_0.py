# Initial grid setup
grid = [
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
]

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Button A controls
button_a_controls = [(1, 2), (1, 0)]

# Button B controls
button_b_controls = [(1, 0), (1, 1), (0, 0)]

# Round 1: Button A
for x, y in button_a_controls:
    toggle_light(grid, x, y)

# Round 2: Button B
for x, y in button_b_controls:
    toggle_light(grid, x, y)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(3) for j in range(3)]

# Print the result
print(result)