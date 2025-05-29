# Initial grid configuration
grid = [
    [1, 1, 0, 1],
    [0, 0, 1, 1],
    [1, 0, 0, 0],
    [1, 1, 1, 1]
]

# Button controls
button_A_controls = [(0, 0), (1, 3), (2, 1), (3, 2), (0, 1), (3, 1), (1, 1)]
button_B_controls = [(0, 2), (3, 1), (0, 1), (2, 1), (3, 2), (1, 1), (2, 2)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'A', 'B']
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A_controls)
    elif round == 'B':
        toggle_lights(grid, button_B_controls)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)