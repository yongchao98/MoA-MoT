# Initial grid configuration
grid = [
    [0, 1, 1],
    [0, 1, 0],
    [1, 0, 1]
]

# Button control mappings
button_A_controls = [(0, 0), (2, 1)]  # 1-based to 0-based index
button_B_controls = [(0, 0), (1, 0)]

# Function to toggle the lights
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

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)