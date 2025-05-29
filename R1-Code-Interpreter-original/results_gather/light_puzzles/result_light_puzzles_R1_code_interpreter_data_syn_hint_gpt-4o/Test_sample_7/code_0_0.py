# Initial grid configuration
grid = [
    [0, 1, 0, 1, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 1],
    [1, 0, 0, 1, 1],
    [0, 1, 1, 1, 1]
]

# Button control mappings
button_A = [(0, 2), (3, 3), (4, 4)]
button_B = [(4, 2), (4, 1), (2, 4), (4, 4), (0, 4), (0, 0)]
button_C = [(0, 0), (3, 1), (2, 2), (3, 3), (4, 2), (3, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

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