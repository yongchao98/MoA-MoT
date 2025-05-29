# Initial grid configuration
grid = [
    [1, 0, 1, 0, 0],
    [1, 1, 0, 1, 1],
    [0, 1, 1, 1, 0],
    [0, 0, 1, 1, 1],
    [1, 1, 1, 1, 1]
]

# Button control mappings
button_A = [(5, 5), (1, 5), (1, 4), (3, 5), (5, 2), (2, 4), (3, 2), (4, 4), (2, 5), (5, 3), (2, 2), (1, 1)]
button_B = [(5, 4), (4, 5), (5, 3)]
button_C = [(4, 1), (4, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = ['C', 'A', 'A', 'A']
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A)
    elif round == 'B':
        toggle_lights(grid, button_B)
    elif round == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)