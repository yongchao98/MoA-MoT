# Initial grid state
grid = [
    [0, 1, 1, 0],
    [1, 1, 1, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 0]
]

# Button control positions
button_A = [(2, 2), (2, 3), (3, 1), (2, 1), (0, 3)]
button_B = [(0, 3), (3, 0), (1, 1), (1, 0)]
button_C = [(3, 0), (3, 2), (0, 0), (1, 0), (0, 1), (0, 2), (1, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_B)  # Round 3: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)