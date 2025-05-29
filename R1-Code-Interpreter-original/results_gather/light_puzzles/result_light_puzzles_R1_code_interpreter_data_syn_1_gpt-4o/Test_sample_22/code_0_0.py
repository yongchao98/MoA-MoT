# Initial grid configuration
grid = [
    [0, 1, 0, 1, 1],
    [1, 0, 1, 0, 0],
    [1, 1, 1, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 0, 1, 1, 1]
]

# Button effects
button_A = [(0, 2), (1, 2)]
button_B = [(3, 2), (3, 1)]
button_C = [(4, 3), (1, 1), (1, 2), (2, 1), (0, 4), (1, 4)]

# Sequence of button presses
button_sequence = ['A', 'B', 'A', 'A']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

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