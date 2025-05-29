# Initial grid configuration
grid = [
    [0, 0, 1, 1],
    [0, 0, 1, 1],
    [1, 0, 1, 0],
    [0, 0, 0, 1]
]

# Button controls
button_A_controls = [(0, 2), (3, 0), (3, 3), (1, 3), (0, 1), (2, 3), (2, 0), (2, 2)]
button_B_controls = [(1, 3), (1, 1)]

# Rounds of button presses
rounds = ['B', 'B', 'B']

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A_controls)
    elif round == 'B':
        toggle_lights(grid, button_B_controls)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)