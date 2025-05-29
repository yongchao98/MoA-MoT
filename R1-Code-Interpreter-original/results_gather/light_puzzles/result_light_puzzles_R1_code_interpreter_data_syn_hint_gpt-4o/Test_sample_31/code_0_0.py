# Initial grid configuration
grid = [
    [0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 1, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 0]
]

# Button control mappings
button_A = [(3, 2), (2, 2), (1, 2), (1, 4), (4, 1)]
button_B = [(2, 5), (2, 4), (4, 5), (2, 1), (2, 3)]
button_C = [(3, 3), (5, 3), (1, 2), (3, 4), (1, 5), (3, 2), (4, 3), (2, 4), (5, 1)]

# Sequence of button presses
button_sequence = ['B', 'C', 'B', 'A']

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the button presses
for button in button_sequence:
    if button == 'A':
        toggle_lights(grid, button_A)
    elif button == 'B':
        toggle_lights(grid, button_B)
    elif button == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)