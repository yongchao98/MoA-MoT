# Initial grid configuration
grid = [
    [1, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 1, 0, 1],
    [0, 0, 1, 1]
]

# Button control mappings
button_A = [(2, 3), (1, 2), (3, 3), (1, 3), (3, 0)]
button_B = [(0, 3), (1, 2), (3, 2), (0, 1), (1, 1), (3, 3)]
button_C = [(0, 1), (3, 3), (3, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_C)  # Round 2: Button C
toggle_lights(grid, button_A)  # Round 3: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)