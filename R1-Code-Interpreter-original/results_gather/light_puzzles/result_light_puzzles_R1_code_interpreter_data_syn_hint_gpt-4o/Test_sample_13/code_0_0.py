# Initial grid configuration
grid = [
    [1, 0, 0, 0, 1],
    [1, 1, 1, 1, 0],
    [1, 1, 0, 1, 1],
    [1, 0, 1, 1, 0],
    [1, 0, 0, 0, 0]
]

# Button control mappings
button_A = [(3, 2), (2, 2), (4, 4), (4, 1), (1, 4), (3, 4), (5, 1), (2, 4), (3, 5)]
button_B = [(3, 4), (3, 5), (2, 5)]
button_C = [(1, 1), (4, 4), (1, 5), (2, 5), (3, 4), (5, 3), (5, 5), (1, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_C)  # Round 1: Button C
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_A)  # Round 3: Button A
toggle_lights(grid, button_B)  # Round 4: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)