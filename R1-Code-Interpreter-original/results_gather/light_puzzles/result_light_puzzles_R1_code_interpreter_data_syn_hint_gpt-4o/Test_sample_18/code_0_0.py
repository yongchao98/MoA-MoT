# Initial grid state
grid = [
    [0, 0, 1, 0, 0],
    [0, 0, 1, 0, 1],
    [1, 1, 0, 1, 0],
    [1, 1, 1, 1, 0],
    [0, 0, 1, 1, 0]
]

# Button control positions
button_A = [(3, 3), (3, 1), (1, 2)]
button_B = [(1, 4), (5, 4), (5, 2), (5, 5), (1, 3), (1, 1), (4, 4), (1, 5), (4, 1), (1, 2), (2, 1), (5, 1)]
button_C = [(1, 2), (4, 5), (3, 2), (3, 1), (1, 3), (2, 3), (3, 5), (2, 5), (4, 1), (5, 3), (1, 5), (3, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_C)  # Round 1: Button C
toggle_lights(grid, button_C)  # Round 2: Button C
toggle_lights(grid, button_B)  # Round 3: Button B
toggle_lights(grid, button_A)  # Round 4: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)