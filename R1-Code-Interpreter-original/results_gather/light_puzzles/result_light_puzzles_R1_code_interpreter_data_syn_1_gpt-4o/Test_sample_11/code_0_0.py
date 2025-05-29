# Initial grid state
grid = [
    [1, 0, 1, 1, 1],
    [0, 1, 0, 0, 1],
    [0, 1, 1, 0, 0],
    [1, 0, 1, 1, 0],
    [1, 1, 1, 0, 1]
]

# Button control positions
button_A = [(3, 5), (2, 3), (2, 2), (4, 5), (4, 2), (1, 4), (3, 4), (1, 1), (3, 3), (2, 5)]
button_B = [(5, 1), (3, 4), (2, 5), (1, 1), (5, 4), (2, 2), (4, 5), (3, 5)]
button_C = [(2, 1), (5, 2), (1, 1), (4, 1), (2, 5), (3, 3), (2, 2), (3, 2), (5, 4), (5, 3), (4, 2), (3, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Sequence of button presses
rounds = [button_C, button_C, button_A, button_B]

# Apply each round
for button in rounds:
    toggle_lights(grid, button)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)