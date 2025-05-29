# Initial grid configuration
grid = [
    [0, 1, 0, 0, 0],
    [1, 1, 1, 0, 0],
    [1, 1, 0, 1, 1],
    [0, 1, 0, 0, 1],
    [1, 1, 1, 1, 0]
]

# Button control mappings
button_A = [(2, 3), (1, 1), (3, 2), (3, 0), (1, 2)]
button_B = [(2, 3), (3, 2), (0, 1), (1, 2), (1, 0), (0, 4), (4, 3), (1, 3), (3, 4), (2, 0), (4, 4)]
button_C = [(0, 3), (3, 2), (2, 0)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_A)  # Round 3: Button A
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)