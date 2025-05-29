# Initial grid configuration
grid = [
    [0, 1, 1],
    [0, 1, 0],
    [1, 0, 1]
]

# Button controls
button_A_controls = [(0, 0), (2, 1)]  # 1-based to 0-based index
button_B_controls = [(0, 0), (1, 0)]

# Function to toggle the lights controlled by a button
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'B']

for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A_controls)
    elif round == 'B':
        toggle_lights(grid, button_B_controls)

# Flatten the grid to the required output format
final_state = [light for row in grid for light in row]
print(final_state)