# Initial grid configuration
grid = [
    [1, 0, 0, 1, 1],
    [0, 0, 1, 1, 0],
    [1, 0, 1, 0, 0],
    [0, 0, 1, 1, 0],
    [1, 0, 1, 1, 0]
]

# Button control mappings
button_A = [(1, 2), (0, 3), (2, 1)]
button_B = [(4, 3), (3, 1), (2, 2), (4, 4), (4, 1), (1, 2), (0, 2), (1, 3), (2, 1), (2, 4), (2, 3)]
button_C = [(2, 3), (4, 4), (4, 3)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_B)  # Round 4: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)