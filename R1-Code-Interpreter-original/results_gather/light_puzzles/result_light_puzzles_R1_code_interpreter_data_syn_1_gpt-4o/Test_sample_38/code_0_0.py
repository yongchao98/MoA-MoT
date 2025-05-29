# Initial grid state
grid = [
    [1, 1, 1, 1, 1],
    [0, 1, 1, 0, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 0],
    [0, 1, 1, 0, 0]
]

# Button control positions
button_A = [(5, 5), (2, 1), (5, 2)]
button_B = [(2, 1), (3, 5), (4, 4), (3, 1), (2, 5), (5, 5), (1, 4), (2, 4)]
button_C = [(2, 2), (3, 1), (2, 5), (2, 1), (3, 5), (4, 5), (2, 4), (4, 4), (5, 3), (1, 4)]

# Convert 1-based positions to 0-based for Python indexing
def convert_positions(positions):
    return [(x-1, y-1) for x, y in positions]

button_A = convert_positions(button_A)
button_B = convert_positions(button_B)
button_C = convert_positions(button_C)

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

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