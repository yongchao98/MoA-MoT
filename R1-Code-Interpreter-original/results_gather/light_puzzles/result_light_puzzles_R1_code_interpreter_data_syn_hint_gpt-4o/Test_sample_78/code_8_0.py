# Initial grid setup
grid = [
    [0, 0, 1, 1],
    [0, 0, 1, 1],
    [1, 0, 1, 0],
    [0, 0, 0, 1]
]

# Button control mappings
button_B_controls = [(1, 3), (1, 1)]

# Function to toggle the lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
for _ in range(3):  # Three rounds of pressing Button B
    toggle_lights(grid, button_B_controls)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)