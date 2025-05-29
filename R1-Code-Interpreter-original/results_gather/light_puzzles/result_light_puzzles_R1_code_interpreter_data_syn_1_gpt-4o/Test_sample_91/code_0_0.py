# Initial grid state
grid = [
    [0, 1, 0, 0],
    [1, 0, 0, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 0]
]

# Button control positions (1-based index)
button_A_controls = [(1, 4), (1, 3), (1, 2), (3, 2)]
button_B_controls = [(1, 4), (3, 1), (2, 3), (4, 1), (2, 1)]

# Convert 1-based index to 0-based index for easier access in Python
button_A_controls = [(x-1, y-1) for x, y in button_A_controls]
button_B_controls = [(x-1, y-1) for x, y in button_B_controls]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_B_controls)  # Round 1: Button B
toggle_lights(grid, button_A_controls)  # Round 2: Button A
toggle_lights(grid, button_A_controls)  # Round 3: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)