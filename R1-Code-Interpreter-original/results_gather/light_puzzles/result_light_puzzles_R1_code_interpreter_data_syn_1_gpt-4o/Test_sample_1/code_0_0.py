# Initial grid state
grid = [
    [1, 0, 0, 0, 0],
    [1, 0, 0, 1, 0],
    [0, 0, 1, 0, 0],
    [1, 1, 0, 0, 1],
    [0, 1, 0, 1, 0]
]

# Button control positions
button_A = [(2, 2), (5, 1), (3, 3), (5, 3), (1, 2), (2, 4), (2, 3)]
button_B = [(2, 2), (1, 2), (5, 4), (4, 2), (4, 1), (3, 2), (2, 3)]
button_C = [(3, 1), (4, 3), (1, 1), (4, 1), (5, 4), (2, 5), (3, 3), (4, 2), (5, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = ['C', 'C', 'C', 'A']
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