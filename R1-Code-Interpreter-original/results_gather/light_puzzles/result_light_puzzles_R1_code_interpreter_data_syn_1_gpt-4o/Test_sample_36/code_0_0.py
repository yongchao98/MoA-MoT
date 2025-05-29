# Initial grid configuration
grid = [
    [0, 1, 1, 0, 0],
    [1, 0, 1, 0, 1],
    [0, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [0, 0, 0, 0, 1]
]

# Button control mappings
button_controls = {
    'A': [(3, 2), (4, 4)],
    'B': [(0, 3), (4, 3), (1, 3), (1, 0), (3, 3)],
    'C': [(2, 1), (4, 1), (4, 0), (0, 4), (1, 2), (0, 2), (4, 2), (3, 2), (1, 4), (1, 0)]
}

# Rounds of button presses
rounds = ['A', 'A', 'B', 'C']

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate each round
for button in rounds:
    toggle_lights(grid, button_controls[button])

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)