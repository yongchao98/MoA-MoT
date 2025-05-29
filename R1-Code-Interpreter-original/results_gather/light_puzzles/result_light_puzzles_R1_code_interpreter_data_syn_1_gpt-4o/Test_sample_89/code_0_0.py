# Initial grid state
grid = [
    [0, 1, 0, 0],
    [1, 0, 1, 1],
    [0, 1, 0, 0],
    [1, 1, 0, 1]
]

# Button control positions (0-indexed)
button_A_positions = [(0, 2), (0, 1), (1, 1), (1, 3), (2, 1), (0, 0)]
button_B_positions = [(0, 1), (2, 3), (2, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A_positions)  # Round 1: Button A
toggle_lights(grid, button_B_positions)  # Round 2: Button B
toggle_lights(grid, button_B_positions)  # Round 3: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)