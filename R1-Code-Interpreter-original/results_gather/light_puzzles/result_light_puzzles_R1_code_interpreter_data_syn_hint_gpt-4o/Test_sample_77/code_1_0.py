# Initial grid configuration
grid = [
    [0, 1, 1],
    [0, 0, 0],
    [0, 0, 1]
]

# Button control mappings
button_A_controls = [(2, 0), (2, 1)]  # 0-indexed positions
button_B_controls = [(2, 1), (2, 0)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A_controls)  # Round 1: Button A
toggle_lights(grid, button_B_controls)  # Round 2: Button B

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)