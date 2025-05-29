# Initial grid state
grid = [
    [1, 1, 1, 1, 0],
    [1, 1, 1, 0, 1],
    [0, 1, 0, 1, 1],
    [1, 0, 1, 0, 0],
    [1, 1, 0, 0, 1]
]

# Button control positions
button_A = [(1, 3), (5, 4), (3, 3), (5, 3), (2, 1)]
button_B = [(1, 1), (2, 4), (2, 2), (1, 4), (5, 3), (5, 5), (4, 5), (3, 3), (3, 5), (5, 2), (3, 2)]
button_C = [(1, 4), (1, 5), (2, 4), (5, 2), (2, 3), (3, 4), (4, 1), (3, 5), (2, 1), (4, 4), (1, 3), (4, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
rounds = [button_B, button_C, button_A, button_A]

for button in rounds:
    toggle_lights(grid, button)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)