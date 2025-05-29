# Initial grid state
grid = [
    [1, 1, 0, 0, 0],
    [1, 1, 0, 0, 0],
    [1, 1, 1, 1, 0],
    [0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0]
]

# Button control positions
button_A = [(5, 4), (5, 1), (5, 2), (1, 3), (2, 3), (3, 1), (1, 2), (5, 3), (3, 2), (1, 4), (4, 2), (5, 5)]
button_B = [(5, 5), (1, 1), (1, 4), (2, 5), (1, 3), (2, 1), (2, 4), (4, 5), (5, 4)]
button_C = [(3, 2), (4, 4), (1, 5), (2, 2), (2, 5), (5, 4)]

# Rounds of button presses
rounds = ['B', 'A', 'B', 'A']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

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