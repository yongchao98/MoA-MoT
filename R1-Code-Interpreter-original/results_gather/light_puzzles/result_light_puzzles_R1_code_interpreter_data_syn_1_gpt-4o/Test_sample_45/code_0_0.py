# Initial grid state
grid = [
    [1, 1, 1, 0, 0],
    [0, 0, 1, 1, 0],
    [1, 0, 0, 0, 1],
    [1, 1, 0, 1, 1],
    [1, 0, 1, 1, 1]
]

# Button control positions
button_A = [(0, 1), (1, 2), (3, 1)]
button_B = [(1, 1), (0, 3), (3, 2), (0, 4), (2, 0), (4, 1), (3, 4), (0, 0), (3, 3), (2, 3)]
button_C = [(1, 2), (0, 2), (1, 4), (2, 3), (0, 0)]

# Rounds of button presses
rounds = ['B', 'A', 'C', 'A']

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