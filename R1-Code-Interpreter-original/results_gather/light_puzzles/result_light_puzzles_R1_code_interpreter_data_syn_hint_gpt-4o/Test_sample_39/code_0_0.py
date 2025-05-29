# Initial grid state
grid = [
    [0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [1, 1, 0, 0, 1],
    [1, 0, 1, 0, 1]
]

# Button control positions (0-indexed)
button_A = [(3, 0), (1, 1)]
button_B = [(2, 2), (4, 1), (0, 4), (4, 2), (1, 2), (0, 0)]
button_C = [(4, 1), (2, 3), (3, 0), (4, 3), (0, 1), (0, 3), (4, 4), (0, 4), (4, 2), (0, 0)]

# Rounds of button presses
rounds = ['C', 'C', 'B', 'B']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate each round
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