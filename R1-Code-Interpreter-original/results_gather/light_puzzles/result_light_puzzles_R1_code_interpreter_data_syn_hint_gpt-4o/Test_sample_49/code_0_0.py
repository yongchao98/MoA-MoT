# Initial grid state
grid = [
    [1, 0, 1, 0, 0],
    [1, 1, 1, 1, 1],
    [1, 0, 1, 1, 0],
    [1, 0, 1, 0, 1],
    [1, 1, 0, 1, 0]
]

# Button control positions (1-based index)
button_A = [(2, 3), (3, 1), (4, 3)]
button_B = [(3, 2), (1, 5), (5, 5), (2, 3), (3, 5), (2, 4), (1, 4), (2, 5), (5, 2), (5, 3), (2, 2)]
button_C = [(2, 2), (1, 4), (4, 1), (5, 4), (1, 3), (3, 4), (1, 2), (5, 1), (5, 2), (2, 3), (2, 5), (3, 3)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'B', 'B', 'A']

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