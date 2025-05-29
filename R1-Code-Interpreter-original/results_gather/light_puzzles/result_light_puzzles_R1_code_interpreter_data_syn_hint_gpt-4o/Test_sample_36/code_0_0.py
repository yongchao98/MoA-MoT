# Initial grid configuration
grid = [
    [0, 1, 1, 0, 0],
    [1, 0, 1, 0, 1],
    [0, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [0, 0, 0, 0, 1]
]

# Button control mappings
button_A = [(3, 2), (4, 4)]
button_B = [(0, 3), (4, 3), (1, 3), (1, 0), (3, 3)]
button_C = [(2, 1), (4, 1), (4, 0), (0, 4), (1, 2), (0, 2), (4, 2), (3, 2), (1, 4), (1, 0)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_B)  # Round 3: Button B
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)