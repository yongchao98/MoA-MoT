# Initial grid setup
grid = [
    [0, 0, 1, 0],
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [1, 0, 1, 1]
]

# Button B controls these positions
button_b_controls = [(1, 0), (1, 2), (2, 0), (2, 3), (1, 3), (2, 2), (2, 1), (3, 2)]

# Function to toggle the lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate pressing Button B three times
for _ in range(3):
    toggle_lights(grid, button_b_controls)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(4) for j in range(4)]

# Print the result
print(result)