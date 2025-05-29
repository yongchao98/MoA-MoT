# Initial grid configuration
grid = [
    [0, 1, 1, 0],
    [0, 1, 1, 1],
    [0, 1, 1, 0],
    [0, 1, 0, 1]
]

# Button A controls these positions
button_A_positions = [(0, 1), (2, 0), (3, 0), (0, 2)]

# Button B controls these positions
button_B_positions = [(3, 1), (1, 1)]

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
# Round 1: Button A
toggle_lights(grid, button_A_positions)

# Round 2: Button B
toggle_lights(grid, button_B_positions)

# Round 3: Button A
toggle_lights(grid, button_A_positions)

# Flatten the grid to a single list for the final output
final_state = [light for row in grid for light in row]
print(final_state)