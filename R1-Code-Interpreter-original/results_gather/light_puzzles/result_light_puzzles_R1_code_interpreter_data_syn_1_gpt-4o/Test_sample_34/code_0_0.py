# Initial grid state
grid = [
    [0, 0, 1, 0, 0],
    [1, 1, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 1, 1, 0]
]

# Button control mappings
button_A = [(3, 1), (3, 0)]
button_B = [(2, 3), (3, 1), (4, 2), (1, 4)]
button_C = [(2, 0), (4, 0), (4, 3)]

# Rounds of button presses
rounds = ['C', 'C', 'B', 'B']

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

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