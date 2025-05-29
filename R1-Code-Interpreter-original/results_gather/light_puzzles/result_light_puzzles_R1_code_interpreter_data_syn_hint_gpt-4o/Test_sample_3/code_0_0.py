# Initial grid state
grid = [
    [1, 1, 1, 0, 1],
    [0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0],
    [0, 1, 1, 0, 1],
    [0, 1, 1, 1, 0]
]

# Button control mappings
button_A = [(1, 3), (2, 4)]
button_B = [(3, 1), (1, 5), (5, 4)]
button_C = [(1, 4), (4, 4), (3, 5), (2, 2), (5, 2), (4, 5), (3, 2), (3, 1), (5, 1), (3, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)