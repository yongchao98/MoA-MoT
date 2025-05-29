# Initial grid configuration
grid = [
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 1],
    [0, 0, 0, 1, 1],
    [1, 1, 1, 0, 0],
    [1, 0, 1, 0, 0]
]

# Button controls
button_A = [(2, 4), (4, 2), (2, 3), (3, 5), (4, 4), (1, 3), (1, 1), (4, 3), (3, 1), (5, 4), (4, 1)]
button_B = [(1, 2), (3, 5), (4, 2), (5, 3), (3, 4), (3, 1), (3, 3)]
button_C = [(5, 5), (5, 1), (4, 2), (4, 1), (3, 1), (2, 3), (3, 4), (5, 2), (1, 5), (1, 1), (3, 5)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = ['C', 'A', 'C', 'C']
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A)
    elif round == 'B':
        toggle_lights(grid, button_B)
    elif round == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]
print(final_state)