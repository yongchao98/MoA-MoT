# Initial grid setup
grid = [
    [0, 0, 1, 0],
    [0, 0, 0, 0],
    [1, 1, 0, 1],
    [0, 0, 0, 1]
]

# Button B controls these positions
button_b_controls = [(0, 1), (0, 0), (1, 0), (3, 1), (1, 1)]

# Function to toggle the lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate pressing Button B three times
for _ in range(3):
    toggle_lights(grid, button_b_controls)

# Flatten the grid to a single list
result = [light for row in grid for light in row]

# Print the result
print(result)