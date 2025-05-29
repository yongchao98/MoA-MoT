# Initial grid configuration
grid = [
    [1, 1, 1, 1, 0],
    [1, 1, 1, 0, 1],
    [0, 1, 0, 1, 1],
    [1, 0, 1, 0, 0],
    [1, 1, 0, 0, 1]
]

# Button controls
button_A = [(1, 3), (5, 4), (3, 3), (5, 3), (2, 1)]
button_B = [(1, 1), (2, 4), (2, 2), (1, 4), (5, 3), (5, 5), (4, 5), (3, 3), (3, 5), (5, 2), (3, 2)]
button_C = [(1, 4), (1, 5), (2, 4), (5, 2), (2, 3), (3, 4), (4, 1), (3, 5), (2, 1), (4, 4), (1, 3), (4, 5)]

# Sequence of button presses
button_sequence = ['B', 'C', 'A', 'A']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

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