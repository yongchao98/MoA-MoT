# Initial grid configuration
grid = [
    [1, 1, 1],
    [1, 1, 1],
    [0, 1, 1]
]

# Button control mappings
button_A_controls = [(1, 1), (1, 2)]
button_B_controls = [(1, 2), (3, 1)]

# Function to toggle the lights controlled by a button
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
# Round 1: Button A
toggle_lights(grid, button_A_controls)

# Round 2: Button A
toggle_lights(grid, button_A_controls)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)