# Initial grid state
grid = [
    [1, 1, 1, 0, 1],
    [0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0],
    [0, 1, 1, 0, 1],
    [0, 1, 1, 1, 0]
]

# Button control positions
button_A = [(0, 2), (1, 3)]
button_B = [(2, 0), (0, 4), (4, 3)]
button_C = [(0, 3), (3, 3), (2, 4), (1, 1), (4, 1), (3, 4), (2, 1), (2, 0), (4, 0), (2, 3)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)