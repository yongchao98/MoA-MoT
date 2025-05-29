# Initial grid state
grid = [
    [0, 0, 1, 1, 0],
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 0],
    [1, 1, 1, 0, 1],
    [1, 0, 1, 0, 0]
]

# Button control positions
button_A = [(2, 1), (2, 4), (0, 2), (1, 4), (3, 4)]
button_B = [(1, 2), (1, 0), (3, 3), (0, 2), (0, 4), (4, 4), (2, 0)]
button_C = [(2, 2), (0, 2), (2, 3), (3, 0), (2, 4)]

# Rounds of button presses
rounds = ['A', 'C', 'C', 'A']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
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