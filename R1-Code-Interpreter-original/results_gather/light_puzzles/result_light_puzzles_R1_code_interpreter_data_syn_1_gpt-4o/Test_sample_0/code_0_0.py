# Initial grid configuration
grid = [
    [0, 1, 0, 0, 1],
    [1, 1, 0, 0, 0],
    [1, 0, 0, 0, 1],
    [0, 1, 0, 1, 1],
    [0, 0, 0, 0, 1]
]

# Button controls
button_A = [(3, 4), (3, 2), (4, 4), (2, 4), (5, 5), (1, 4), (5, 2), (4, 2), (4, 1), (4, 3)]
button_B = [(1, 3), (1, 5), (2, 3), (1, 2), (1, 4), (4, 4), (2, 1), (3, 2), (5, 3)]
button_C = [(1, 1), (2, 3)]

# Convert 1-based to 0-based indices
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
rounds = [button_B, button_C, button_A, button_C]

for button in rounds:
    toggle_lights(grid, button)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)