# Initial grid state
grid = [
    [1, 1, 0, 0, 1],
    [1, 0, 1, 1, 1],
    [0, 1, 1, 0, 0],
    [1, 1, 0, 1, 0],
    [1, 0, 1, 1, 1]
]

# Button control mappings
button_A = [(5, 2), (4, 4)]
button_B = [(1, 3), (2, 1), (1, 4)]
button_C = [(5, 2), (4, 3), (1, 3), (1, 1), (2, 1), (3, 5), (3, 4), (4, 2), (3, 2), (3, 1), (2, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_A)  # Round 4: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)