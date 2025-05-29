# Initial grid configuration
grid = [
    [1, 0, 1, 1],
    [0, 1, 0, 1],
    [0, 0, 0, 1],
    [0, 1, 1, 1]
]

# Button effects
button_A = [(1, 1), (3, 1)]  # 0-indexed positions
button_B = [(1, 1), (3, 1), (0, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate the rounds
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_B)  # Round 3: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)