# Initial grid state
grid = [
    [0, 0, 1, 0, 1],
    [0, 1, 1, 0, 1],
    [0, 0, 1, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 0, 0, 0, 1]
]

# Button control positions (1-based index)
button_A = [(4, 4), (1, 4), (4, 3), (5, 4), (5, 1), (1, 5)]
button_B = [(3, 3), (5, 2), (2, 3), (4, 2), (4, 5), (5, 1), (5, 4), (3, 2), (2, 5), (3, 5), (3, 1), (2, 1)]
button_C = [(4, 4), (1, 3), (5, 1), (3, 3), (2, 1), (3, 5)]

# Convert 1-based index to 0-based index
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
rounds = ['A', 'C', 'C', 'A']
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