# Initial grid configuration
grid = [
    [0, 1, 1, 0],
    [0, 1, 1, 1],
    [0, 1, 1, 1],
    [0, 0, 0, 0]
]

# Button control mappings
button_A = [(3, 1), (1, 3), (2, 0), (1, 1), (0, 2), (0, 1), (2, 3), (1, 0)]
button_B = [(1, 0), (2, 0), (0, 0)]
button_C = [(0, 1), (3, 3), (2, 3), (2, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the button presses
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)