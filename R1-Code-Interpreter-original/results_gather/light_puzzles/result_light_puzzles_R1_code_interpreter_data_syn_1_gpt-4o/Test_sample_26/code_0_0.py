# Initial grid state
grid = [
    [0, 1, 1, 1, 0],
    [1, 0, 1, 0, 0],
    [1, 1, 0, 1, 1],
    [0, 1, 0, 1, 1],
    [0, 0, 1, 1, 1]
]

# Button control positions
button_A = [(2, 1), (4, 4), (4, 5), (3, 5), (5, 1), (1, 3), (3, 4)]
button_B = [(5, 2), (5, 3), (3, 5), (1, 2), (4, 5), (2, 1), (1, 1), (3, 2), (2, 3), (5, 1)]
button_C = [(3, 4), (1, 4), (2, 2), (4, 4), (1, 2), (2, 5), (2, 4), (4, 1), (3, 1), (2, 1), (1, 5), (3, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
rounds = ['C', 'A', 'A', 'C']
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