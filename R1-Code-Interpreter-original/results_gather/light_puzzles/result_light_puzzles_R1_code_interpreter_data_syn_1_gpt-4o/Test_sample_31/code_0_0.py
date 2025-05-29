# Initial grid configuration
grid = [
    [0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 1, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 0]
]

# Button controls
button_controls = {
    'A': [(3, 2), (2, 2), (1, 2), (1, 4), (4, 1)],
    'B': [(2, 5), (2, 4), (4, 5), (2, 1), (2, 3)],
    'C': [(3, 3), (5, 3), (1, 2), (3, 4), (1, 5), (3, 2), (4, 3), (2, 4), (5, 1)]
}

# Rounds of button presses
rounds = ['B', 'C', 'B', 'A']

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
for button in rounds:
    toggle_lights(grid, button_controls[button])

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)