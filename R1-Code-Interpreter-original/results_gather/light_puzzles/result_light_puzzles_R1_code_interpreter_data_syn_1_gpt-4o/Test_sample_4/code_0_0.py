# Initial grid state
grid = [
    [0, 1, 1, 1, 0],
    [1, 0, 0, 0, 1],
    [0, 0, 0, 1, 0],
    [1, 1, 1, 1, 1],
    [0, 0, 1, 0, 0]
]

# Button control positions (1-based index)
button_A_positions = [(4, 5), (1, 3), (5, 2), (4, 1), (2, 4), (2, 2), (4, 2), (2, 1), (4, 3), (3, 1)]
button_B_positions = [(1, 2), (4, 2), (3, 1), (3, 2), (3, 3), (4, 1), (5, 1), (3, 5), (2, 5)]
button_C_positions = [(2, 5), (2, 2), (4, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = ['B', 'B', 'B', 'C']
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