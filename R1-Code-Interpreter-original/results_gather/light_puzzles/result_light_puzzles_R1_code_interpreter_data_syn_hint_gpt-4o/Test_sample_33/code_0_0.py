# Initial grid configuration
grid = [
    [1, 1, 0, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 1, 1, 0, 1],
    [0, 0, 0, 0, 1],
    [1, 1, 0, 1, 0]
]

# Button controls
button_A = [(1, 2), (4, 3), (5, 2), (3, 3), (4, 1), (5, 3), (4, 4), (5, 5), (2, 4)]
button_B = [(5, 2), (1, 5), (3, 5), (4, 4), (4, 1), (5, 5)]
button_C = [(2, 4), (2, 3), (5, 5)]

# Sequence of button presses
rounds = ['C', 'A', 'A', 'B']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
for button in rounds:
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