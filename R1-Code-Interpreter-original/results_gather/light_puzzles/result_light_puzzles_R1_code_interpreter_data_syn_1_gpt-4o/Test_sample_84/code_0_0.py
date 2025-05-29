# Initial grid configuration
grid = [
    [0, 1, 0, 0],
    [1, 1, 1, 0],
    [0, 1, 1, 0],
    [0, 1, 0, 0]
]

# Button A controls these positions
button_A_positions = [(1, 3), (0, 1), (1, 2), (2, 1)]

# Button B controls these positions
button_B_positions = [(2, 2), (2, 1), (0, 3), (2, 0)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the button presses
for _ in range(3):  # Three rounds of pressing Button B
    toggle_lights(grid, button_B_positions)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)