# Initial grid state
grid = [
    [1, 1, 0, 1, 1],
    [1, 0, 0, 1, 1],
    [1, 1, 0, 1, 1],
    [0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1]
]

# Button control positions (1-based index)
button_A = [(5, 3), (3, 5), (4, 4), (1, 4), (2, 3), (2, 1), (4, 1), (1, 5), (5, 4), (3, 1), (5, 5), (1, 1)]
button_B = [(5, 5), (3, 3), (1, 1), (3, 4), (1, 4), (2, 1), (3, 1), (2, 2)]
button_C = [(3, 5), (1, 5), (4, 2), (5, 3)]

# Convert 1-based index to 0-based index for easier manipulation
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Sequence of button presses
rounds = ['A', 'B', 'B', 'A']

# Apply each round
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