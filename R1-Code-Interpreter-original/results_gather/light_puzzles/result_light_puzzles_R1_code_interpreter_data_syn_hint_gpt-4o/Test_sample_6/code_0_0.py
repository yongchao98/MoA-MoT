# Initial grid state
grid = [
    [0, 1, 0, 1, 0],
    [1, 0, 1, 1, 1],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 1, 1, 1, 0]
]

# Button control positions (1-based index)
button_A = [(1, 4), (1, 3), (5, 1), (4, 5), (3, 1), (4, 3), (5, 5), (3, 4), (5, 2)]
button_B = [(4, 5), (1, 4), (3, 4), (5, 3), (3, 2), (3, 5), (1, 1), (5, 4), (2, 1), (2, 4), (4, 4)]
button_C = [(1, 4), (3, 4), (1, 3), (3, 1), (5, 5), (1, 2), (3, 3), (5, 1), (4, 1)]

# Convert 1-based index to 0-based index for easier manipulation
def convert_to_zero_based(positions):
    return [(x-1, y-1) for x, y in positions]

button_A = convert_to_zero_based(button_A)
button_B = convert_to_zero_based(button_B)
button_C = convert_to_zero_based(button_C)

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'C', 'B', 'C']
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