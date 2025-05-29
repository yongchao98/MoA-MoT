# Initial grid configuration
grid = [
    [1, 0, 1, 1, 0],
    [1, 1, 0, 0, 1],
    [0, 0, 1, 0, 1],
    [1, 1, 1, 0, 1],
    [1, 1, 0, 0, 1]
]

# Button controls
button_A = [(3, 3), (5, 2), (4, 1), (3, 4), (3, 2), (2, 1), (5, 1)]
button_B = [(2, 2), (1, 4), (4, 1), (4, 2), (5, 2), (3, 4), (1, 2), (1, 1), (2, 5), (4, 4)]
button_C = [(3, 2), (5, 4), (1, 3), (2, 1), (4, 4), (2, 4), (1, 2), (2, 5), (5, 5)]

# Sequence of button presses
sequence = ['C', 'A', 'B', 'C']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the button presses
for button in sequence:
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