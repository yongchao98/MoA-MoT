# Initial grid configuration
grid = [
    [0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [1, 0, 1, 1, 0],
    [1, 0, 1, 1, 0],
    [0, 1, 0, 0, 0]
]

# Button control mappings
button_A = [(5, 4), (2, 2), (1, 4), (4, 2), (3, 2), (5, 1), (4, 4)]
button_B = [(4, 3), (4, 5)]
button_C = [(1, 3), (4, 5), (3, 1), (5, 1), (1, 5), (5, 3), (2, 1), (2, 5), (3, 2)]

# Convert 1-based to 0-based indices
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
rounds = [button_A, button_C, button_A, button_B]

for button in rounds:
    toggle_lights(grid, button)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)