# Initial grid state
grid = [
    [1, 1, 0, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 1, 1, 0, 1],
    [0, 0, 0, 0, 1],
    [1, 1, 0, 1, 0]
]

# Button control positions (1-based index)
button_A_positions = [(1, 2), (4, 3), (5, 2), (3, 3), (4, 1), (5, 3), (4, 4), (5, 5), (2, 4)]
button_B_positions = [(5, 2), (1, 5), (3, 5), (4, 4), (4, 1), (5, 5)]
button_C_positions = [(2, 4), (2, 3), (5, 5)]

# Convert 1-based index to 0-based index for easier access in Python
def convert_positions(positions):
    return [(x-1, y-1) for x, y in positions]

button_A_positions = convert_positions(button_A_positions)
button_B_positions = convert_positions(button_B_positions)
button_C_positions = convert_positions(button_C_positions)

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Sequence of button presses
rounds = ['C', 'A', 'A', 'B']

# Apply each round
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A_positions)
    elif round == 'B':
        toggle_lights(grid, button_B_positions)
    elif round == 'C':
        toggle_lights(grid, button_C_positions)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)