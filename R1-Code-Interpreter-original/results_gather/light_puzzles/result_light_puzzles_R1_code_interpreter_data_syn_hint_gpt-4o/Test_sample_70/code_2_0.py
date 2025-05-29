# Initial grid configuration
grid = [
    [0, 0, 0],
    [0, 1, 0],
    [1, 0, 1]
]

# Button controls
button_A_controls = [(1, 0), (1, 1), (2, 1)]
button_B_controls = [(1, 1), (2, 0), (0, 1)]

# Function to toggle lights
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
result = [grid[i][j] for i in range(3) for j in range(3)]
print(result)