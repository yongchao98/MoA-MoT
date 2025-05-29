# Initial grid setup
grid = [
    [1, 0, 1, 0],
    [0, 1, 1, 0],
    [1, 0, 1, 0],
    [1, 1, 0, 1]
]

# Button control positions (1-based index)
button_A_positions = [(4, 4), (3, 4), (2, 2), (1, 3), (2, 4), (1, 4)]
button_B_positions = [(1, 1), (4, 2), (1, 2), (4, 1)]

# Convert 1-based index to 0-based index for easier list manipulation
button_A_positions = [(x-1, y-1) for x, y in button_A_positions]
button_B_positions = [(x-1, y-1) for x, y in button_B_positions]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'B', 'A']

for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A_positions)
    elif round == 'B':
        toggle_lights(grid, button_B_positions)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)