# Initial grid configuration
grid = [
    [0, 0, 1, 0],
    [0, 0, 0, 0],
    [1, 1, 0, 1],
    [0, 0, 0, 1]
]

# Button A controls these positions
button_A_positions = [(1, 0), (3, 0), (0, 3), (1, 1), (0, 0), (1, 3), (2, 2)]

# Button B controls these positions
button_B_positions = [(0, 1), (0, 0), (1, 0), (3, 1), (1, 1)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the button presses
for _ in range(3):  # Round 1, 2, 3: Button B
    toggle_lights(grid, button_B_positions)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)