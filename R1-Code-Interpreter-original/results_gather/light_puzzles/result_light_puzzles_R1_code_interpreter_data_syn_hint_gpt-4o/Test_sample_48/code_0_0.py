# Initial grid configuration
grid = [
    [0, 0, 1, 0, 1],
    [0, 0, 1, 0, 1],
    [1, 1, 0, 0, 1],
    [1, 0, 1, 1, 1],
    [1, 0, 0, 1, 1]
]

# Button controls
button_A = [(1, 5), (5, 3), (2, 3), (5, 2), (4, 2), (5, 5), (1, 1), (1, 4), (4, 3), (1, 3), (2, 1), (2, 2)]
button_B = [(2, 4), (5, 4), (3, 1), (1, 5), (1, 3), (5, 2), (2, 2), (1, 2), (1, 1), (4, 1)]
button_C = [(5, 5), (4, 4), (1, 3), (5, 4), (5, 1), (3, 4), (1, 4), (2, 1), (4, 3), (3, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
rounds = ['C', 'C', 'B', 'B']
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