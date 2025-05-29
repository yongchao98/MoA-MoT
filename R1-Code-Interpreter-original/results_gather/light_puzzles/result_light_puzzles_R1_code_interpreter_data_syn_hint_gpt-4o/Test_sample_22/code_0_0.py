# Initial grid configuration
grid = [
    [0, 1, 0, 1, 1],
    [1, 0, 1, 0, 0],
    [1, 1, 1, 0, 1],
    [1, 0, 1, 0, 0],
    [0, 0, 1, 1, 1]
]

# Button control mappings
button_A = [(0, 2), (1, 2)]
button_B = [(3, 2), (3, 1)]
button_C = [(4, 3), (1, 1), (1, 2), (2, 1), (0, 4), (1, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_A)  # Round 3: Button A
toggle_lights(grid, button_A)  # Round 4: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)