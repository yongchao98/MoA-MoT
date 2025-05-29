# Initial grid configuration
grid = [
    [1, 0, 0, 0],
    [1, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 0]
]

# Button A controls
button_A_controls = [(2, 3), (1, 2)]

# Button B controls
button_B_controls = [(1, 0), (2, 0), (0, 3), (0, 1), (2, 1)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
# Round 1: Button B
toggle_lights(grid, button_B_controls)

# Round 2: Button A
toggle_lights(grid, button_A_controls)

# Round 3: Button A
toggle_lights(grid, button_A_controls)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)